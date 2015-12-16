function [e,data,raw] = xsim(input)                 
%   [e,data,raw] = XSIM(input) solves a stochastic differential equation.
%   It returns averages as defined by the 'input' cell array. 
%   Returned value e is the maximum error vector e = [e(1),e(2)].
%   Here e(1) is the time-step, e(2) the sampling error. 
%   XSIM also calls parallel loops, estimates errors, writes files.
%   The output array 'data' includes the averages and error bars.
%   The output array 'raw' is the raw trajectories if required.
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

tic();                                           %%set timer
input = xinpreferences(input);                   %%set input defaults
latt = xlattice(input);                          %%initialise lattice
r = latt{1};                                     %%common input data{s}
s = RandStream('CombRecursive','Seed',1);        %%define RNG type
RandStream.setGlobalStream(s);                   %%set RNG type
fprintf ('\n \n %s  %s \n',r.version,r.name);    %%integration number
fprintf (' xsim starting\n');                    %%integration start
if r.print >1                                    %%if print switch is verbose
    r                                            %%display lattice data
end;                                             %%end if print switch
sequence = length(input);                        %%find sequence length
for s=1:sequence                                 %%loop over sequence
    data{s} =  zeros(latt{s}.d.data);            %%initialise output data                                        
end                                              %%end loop over sequence
if r.ensembles(3) > 1                            %%if parallel simulations
    parfor npar = 1:r.ensembles(3)               %%loop parallels
        [dp,raw{npar}] = xensemble(npar,latt);   %%call ensembles loop
        data = addcell(data,dp);                 %%call ensembles loop
    end                                          %%end parallels loop
else                                             %%no parallels
    [data,raw] = xensemble(1,latt);              %%call ensembles loop
end                                              %%end if parallels
fprintf('Integration time = %f \n',toc());       %%print time taken
for s=1:sequence                                 %%loop over sequence
  data{s}(3,:,:,:)=data{s}(3,:,:,:) - data{s}(1,:,:,:).^2; 
  e(1)=0;                                        %%initial step errors
  if r.errorchecks >1                            %%if errorchecks
    data{s}(2,:,:,:)=data{s}(2,:,:,:) - data{s}(1,:,:,:);     %%store errors
    if r.order > 0                               %%if extrapolation order
        data{s}(2,:,:,:) = data{s}(2,:,:,:)/(2^(r.order)-1.); %%Error-bar
        data{s}(1,:,:,:) = data{s}(1,:,:,:)- data{s}(2,:,:,:);%%Extrapolate 
    end                                          %%end if extrapolation
    data{s}(2,:,:,:) = abs(data{s}(2,:,:,:));    %%absolute value of error
    e(1) = max(max(max(data{s}(2,:,:,:))));      %%max step error
    fprintf('Max step error = %e\n',e(1));       %%print max step error
  end                                            %%end if errorchecks
  e(2) =0;                                       %%initial sampling errors
  if r.n.ensemble>1                              %%if ensemble averaging
    data{s}(3,:,:,:)=real(sqrt(data{s}(3,:,:,:)/(r.n.ensemble-1)));%%std. dev.
    e(2) = max(max(max(abs(data{s}(3,:,:,:))))); %%max of sd matrix
    fprintf('Max sampling error = %e\n',e(2));   %%sampling error
  end                                            %%end if 
end                                              %%end for seqence                                          
if ~strcmp(r.file,'')                            %%if filename not blank
    xwrite(e,input,data,raw);                    %%save data in file
end                                              %%end if file not blank
end                                              %%end xsim


function data1 = addcell(data1,data2)
sequence = length(data1);
for s = 1:sequence
    data1{s} = data1{s}+data2{s};
end
end

%%Change history: v0.5
%%Include default and calculated values in al struct
%%Print summary data{s}
%%Change history: v0.51
%%Change dimension of data{s}
%%Change history: v0.6
%%Allows for sequenced inputs
%%Change history: v0.61
%%Produces error output for batch testing
%%Change history: v0.7
%%Summaries printed elsewhere, shortened Matlabfile --> Matfile
%%Change history: v0.8
%%Combines data{s} and input in one cell