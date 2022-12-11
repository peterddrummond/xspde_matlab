function [err,data,input,raw] = xsim(input)                 
%   [error,inputdata,raw] = XSIM(input) 
%   Solves stochastic differential and partial differential equations.
%   It returns averages as defined by the 'input' cell array. 
%   XSIM also calls parallel loops, estimates errors, writes files.
%   The 'error' is a six dimensional vector:
%    e(1) is the overall maximum error
%    e(2) is the overall maximum (relative) step error
%    e(3) is the overall maximum (relative) sampling error
%    e(4) is the overall maximum (relative) difference error
%    e(5) is the overall maximum chi-squared error
%    e(6) is the elapsed time
%
%   The 'data' cell array   gives the averages, errors, and comparisons.
%   The 'input' cell array  gives input data including all default values.
%   The 'raw'  cell array   gives raw trajectories if required.
%
%   Called by: xspde
%   Needs: xmakecell, xpreferences, xlattice, xensemble, xaddcell, xwrite
%   Licensed by Peter D. Drummond, (2021) - see License.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% MAKE INPUT A CELL IF NEEDED
%
tic();                                           %%set timer
input = xmakecell(input);                        %%make input a cell
input = xpreferences(input);                     %%get 'input' defaults
sequence = length(input);                        %%get sequence length
data = cell(sequence);                           %%allocate data as cell
l = cell(sequence);                              %%allocate lattice as cell
err = zeros(1,6);                              %%initialize error vector

%%Loop over sequence to initialize data sets

totn =0;
for s=1:sequence                                 %%loop over sequence
  l{s} = xlattice(input{s});                     %%initialise lattice
  if l{s}.seed >=0  && ~l{s}.octave              %%check if not octave
    sd = RandStream('CombRecursive','Seed',1);   %%define RNG type
    RandStream.setGlobalStream(sd);              %%set RNG type
  end 
  for n=1:l{s}.functions
        data{s}{n} = zeros(l{s}.d.data{n});     %%initialise data cells
        totn = totn +1;
  end
  if l{s}.print > 1                              %%if print switch verbose
    fprintf('\n xSIM parameters\n');             %%display input data
    display (l{s});                              %%display input data
  end                                            %%end if print switch
end                                              %%end loop over sequence
p = l{1};                                        %%first sequence cell
xpr(1,p,'\nxSIM starting, %s\n\n',input{1}.version);%%start up message

%%Carry out parallel simulations if needed

if p.ensembles(3) > 1                            %%if parallel simulations
    parfor npar = 1:p.ensembles(3)               %%loop on parallel threads
        tic();                                   %%set timer
        [dp,raw1] = xensemble(npar,l);           %%call xensemble function
        raw(:,:,:,npar) =raw1;                   %%accumulate raw data
        data = xaddcell(data,dp);                %%accumulate averages
    end                                          %%end parallels loop
    raw = reshape(raw,[sequence,p.errorchecks,p.ncopies]); %%reshape raw 
else                                             %%no parallels specified
    tic();                                       %%set timer
    [data,raw] = xensemble(1,l);                 %%call xensemble function
end                                              %%end if parallels

%%Loop over sequence to finalize data sets and error summaries

echisq = 0;                                      %%initial chi-square
compsq = 0;                                      %%Total chi-square points
ns = 0;                                          %%total observe number
stepn = zeros(1,totn);                           %%initialize step errors
diffn = zeros(1,totn);                           %%initialize diff errors
sampn = zeros(1,totn);                           %%initialize sample error
for s=1:sequence                                 %%loop over sequence
  p = l{s};                                      %%get lattice for sequence
  xpr(1,p,'\nDataset %d: %s\n',s,p.name);        %%print name for sequence
  ens = prod(p.ensembles);                       %%initialize ensembles
  p.bin{1+p.functions}=[];                       %%initialize bin arrays
  
%%Loop over functions of data sets and errors

  for n=1:p.functions                            %%loop on graph functions
    ns =  ns+1;                                  %%initialize count                     
    chisum = 0;                                  %%initialize chi-sq errors
    d=data{s}{n};                                %%get n-th data array
    en = max(abs(d(:,1)));                       %%get error normalisation
    rel = 'relative';                            %%relative normalisation
    if p.ncopies > 1                             %%if higher ensembles
        d(:,3) = d(:,3) - abs(d(:,1)).^2;        %%get variance
        d(:,3)=abs(sqrt(d(:,3)/(p.ncopies-1)));  %%get stand. dev. in mean
    else                                         %%no higher ensembles
        d(:,3) = 0.0*d(:,3);                     %%set stand. dev. to zero
    end                                          %%end if ensemble
    if p.bins{n} > 0                             %%If probabilities needed
       p.scale{n} = p.a.bins{n}*ens;             %%scale variance estimates
       d(:,3) = sqrt(d(:,1)/p.scale{n});         %%Poissonian sd estimate
       p.cutoffs{n} = max(p.cutoffs{n},p.mincount/p.scale{n});
    end                                          %%end if probabilities
    sampn(ns) = sqrt(mean(d(:,3)'.^2));          %%n-th rms sampling error 
    if p.errorchecks > 1                         %%if errorchecks needed
      d(:,2) = d(:,2) - d(:,1);                  %%get difference in result
      if p.order > 0                             %%if extrapolation order
        d(:,2) = d(:,2)/(2^(p.order)-1.);        %%estimate true error
        d(:,1) = d(:,1) - d(:,2);                %%extrapolate
      end                                        %%end if extrapolation
      d(:,2) = abs(d(:,2));                      %%get the absolute error
      stepn(ns) = sqrt(mean(d(:,2).^2));         %%n-th rms step error  
    end                                          %%end if errorchecks
    d = reshape(d,p.gpoints{n});                 %%output: full indexing
    for nd=2:p.dimension                         %%loop over dimensions
      if p.ftransforms{n}(nd)                    %%check Fourier transforms
         d = fftshift(d,1+nd);                   %%shift Fourier origin
      end                                        %%end check Fourier
    end                                          %%end loop over dimensions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% GET COMPARISONS AND COMPRESS
%
    [d,p.xk{n},p.axes{n}] = xcompress(n,d,p);    %%Compress data to axes    
    p = p.grid(n,p);                             %%get graphics x,o,grid
    p.compares{n} = 0;                           %%initialize compares
    if ~isempty(p.compare{n})                    %%If comparisons needed
      nd = length(p.gpoints{n});                 %%Get data dimensions
	  c = p.compare{n}(p);
      c = c+zeros(p.gpoints{n}(1:end-1));        %%extend compare function
      szc = size(c);                             %%Size of compare data
      comperrors = 1;                            %%Compare errors switch
      if (length(szc) < nd)                      %%If compare has no errors
         c = cat(nd,c,0*c,0*c);                  %%Extend compare size
         comperrors = 0;                         %%No compare errors
      end                                        %%End if compare
      p.compares{n} = 3;                         %%Reset size of compares
    end                                          %%End if comparisons
    dsize = size(d);                             %%Get compressed data size
    lastd = length(dsize);                       %%Get last data dimension
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PERFORM CHI-SQUARE/DIFF TESTS
%    
    if p.compares{n} > 0                         %%If comparisons needed
      dsize=size(d);                             %%Get data size
      c = xcompress(n,c,p);                      %%Compress data to axes
      csize = size(c);                           %%Comparison size
      cpoint = prod(dsize(1:(lastd-1)));         %%n-th data points
      d = reshape(d,[cpoint,3]);                 %%Flatten data
      c = reshape(c,[cpoint,p.compares{n}]);     %%Get comparisons
      if ~isempty(p.bin{n}) && ~comperrors       %%If probabilities needed
             d(:,3) = sqrt(c(:,1)/p.scale{n});   %%Use analytic sd estimate
      end                                        %%End if probabilities
      diff = abs(d(:,1)-c(:,1));                 %%Comparison differences
      if max(abs(c(:,1))) < 1.e-100              %%If comparison too small
          en = 1;                                %%Don't normalise 
          rel = 'absolute';                      %%absolute normalisation
      end                                        %%End if comparison small         
      sig = d(:,2).^2 + d(:,3).^2;               %%Data error squared
      for ind = 2:p.compares{n}                  %%For comparison errors
          sig = sig + c(:,ind).^2;               %%Comparison error squared
      end                                        %%End for comp. errors
      sigm = max(sig);                           %%Maximum errors
      diffn(ns) = sqrt(mean(diff.^2))/en;        %%RMS differences
      epoint = sig > 0;                          %%Significant error points
      if p.cutoffs{n}                            %%If there are cutoffs
        epoint = epoint&(c(:,1)>p.cutoffs{n});   %%Check for cutoffs in c
        epoint = epoint&(d(:,1)>p.cutoffs{n});   %%Check for cutoffs in d
      end                                        %%End if there are cutoffs
      sigp = sum(epoint);                        %%Count significant points
      if sigp > 0 && sigm > 0 && (p.ncopies > 1 || ~isempty(p.bin{n}))
        chisq = diff.^2./( 1.e-20+sig);          %%Relative error squares
        chisum = sum(chisq.*epoint);             %%Get chi-square  
        echisq = echisq+chisum;                  %%Total chisq errors
        compsq =  compsq+sigp;                   %%Total chisq points
      end                                        %%End if chi-square
      d = reshape(d,dsize);                      %%Unflatten data
      c = reshape(c,csize);                      %%Unflatten comparisons
      d = cat(lastd,d,c);                        %%Concatenate with data
    end                                          %%End if comparisons
    xpr(1,p,'\nFunction %d: %s\n',n,p.olabels{n});%%print function
    if stepn(ns) > 0                             %%if step error
      stepn(ns) = stepn(ns)/en;                  %%n-th RMS step error            
      xpr(1,p,'RMS time-step error  = %.3g (%s)\n',stepn(ns),rel);
    end  
    if sampn(ns) > 0                             %%if sampling error
      sampn(ns) = sampn(ns)/en;                  %%n-th RMS step error 
      xpr(1,p,'RMS sampling error   = %.3g (%s)\n',sampn(ns),rel);
    end
    if diffn(ns) > 0                             %%if diff error   
         xpr(1,p,'RMS comparison error = %.3g (%s)\n',diffn(ns),rel);
    end
    if chisum > 0                 %%if chi-square
      err(5) = err(5) + chisum;
      err(6) = err(6) + sigp;
      xpr(1,p,'Chi-square/points    = %.3g\n',chisum/sigp);
    end
    data{s}{n} = d;                              %%Store data in cell array
  end                                            %%end loop on functions
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% PRINT SUMMARIES AND STORE
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% DELETE TEMPORARY ARRAYS, STORE r
%
  p = rmfield(p,{'r','k','x','y','z','kx','ky','kz','t','w'});
  if p.dimension > 1
    p = rmfield(p,{'D','Dx','Dy','Dz','boundval','infilt','noisefilt'});
  end
  input{s}=p;                                    %%store new parameters
end                                              %%end loop on seqence
err(2) = sqrt(sum(stepn.^2.)/(1.e-99+sum(stepn>0))); %%store error
err(3) = sqrt(sum(sampn.^2.)/(1.e-99+sum(sampn>0))); %%store error
err(4) = sqrt(sum(diffn.^2.)/(1.e-99+sum(diffn>0))); %%store error
%
err(1) = sqrt(sum(err(2:4).^2.)/...
               (1.e-99+sum(err(2:4)>0)));      %%store error
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% WRITE DATA TO FILE AND RETURN
%
if err(5)>0
    err(5) = err(5)/err(6);
end
err(6) = toc(); 
  if err(1)>0
      xpr(0,p,'\nxSIM integration errors for %s\n',input{1}.name);
      xpr(0,p,'RMS overall relative error    = %.3g\n',err(1));
  end
  if err(2)>0    
      xpr(0,p,'RMS relative step error       = %.3g\n',err(2));
  end
  if err(3)>0      
      xpr(0,p,'RMS relative sampling error   = %.3g\n',err(3));
  end 
  if err(4)>0
      xpr(0,p,'RMS relative comparison error = %.3g\n',err(4));
  end
  if err(5)>0
      xpr(0,p,'Chi-square/points             = %.3g\n',err(5));
  end
  xpr(0,p,'\nxSIM simulated %d functions, time taken = %.3gs\n',ns,err(6));
if ~strcmp(p.file,'')                            %%if filename exists
    xwrite(err,data,input,raw);                  %%save data in file
end                                              %%end if file not blank
if length(input) == 1                            %%if  input length is 1
    input = input{1};                            %%returns struct, not cell
end                                              %%end if input length is 1
end                                              %%end xsim
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END XSIM
%
function data1 = xaddcell(data1,data2)
%   data1 = XADDCELL(data1,data2) 
%   Adds two cell arrays together
%   Licensed by Peter D. Drummond, (2015) - see License.txt 

for seq = 1:length(data1)
  for n = 1:length(data1{seq})
    data1{seq}{n} = data1{seq}{n}+data2{seq}{n};
  end
end
end
