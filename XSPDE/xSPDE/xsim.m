function [error,input,data,raw] = xsim(input)                 
%   [error,input,data,raw] = XSIM(input) 
%   Solves stochastic differential and partial differential equations.
%   It returns averages as defined by the 'input' cell array. 
%   Returned value error is the maximum error vector, [error(1),error(2)].
%   Here error(1) is the time-step, error(2) the sampling error. 
%   XSIM also calls parallel loops, estimates errors, writes files.
%   The output cell array 'input' is input data including default values.
%   The output cell array 'data' includes the averages and error bars.
%   The output cell array 'raw' is the raw trajectories if required.
%   MIT license by Peter D. Drummond, (2015) - see License.txt 

tic();                                           %%set timer
input = xinpreferences(input);                   %%set input defaults
latt = xlattice(input);                          %%initialise lattice
r = latt{1};                                     %%first sequence input
s = RandStream('CombRecursive','Seed',1);        %%define RNG type
RandStream.setGlobalStream(s);                   %%set RNG type
fprintf ('\nxSIM%s: %s\n',r.version,r.name);       %%integration name
if r.print >1                                    %%if print switch is verbose
    r                                            %%display lattice data
end;                                             %%end if print switch
sequence = length(input);                        %%find sequence length
for s=1:sequence                                 %%loop over sequence
    data{s} =  zeros(latt{s}.d.data);            %%initialise output data                                        
end                                              %%end loop over sequence
if r.ensembles(3) > 1                            %%if parallel simulations
    parfor npar = 1:r.ensembles(3)               %%loop over parallel threads
        [dp,raw{npar}] = xensemble(npar,latt);   %%call xensemble function
        data = xaddcell(data,dp);                %%accumulate averages
    end                                          %%end parallels loop
else                                             %%no parallels specified
    [data,raw] = xensemble(1,latt);              %%call xensemble function
end                                              %%end if parallels
fprintf('Integration time = %f \n',toc());       %%print time taken
for s=1:sequence                                 %%loop over sequence
  r = latt{s};                                   %%get lattice for this sequence
  data{s}(3,:,:,:)=data{s}(3,:,:,:) - data{s}(1,:,:,:).^2; %%get sampling variance
  error(1)=0;                                    %%initial step errors
  if r.errorchecks > 1                           %%if errorchecks needed
    data{s}(2,:,:,:)=data{s}(2,:,:,:) - data{s}(1,:,:,:);     %%store errors
    if r.order > 0                               %%if extrapolation order
        data{s}(2,:,:,:) = data{s}(2,:,:,:)/(2^(r.order)-1.); %%Error-bar
        data{s}(1,:,:,:) = data{s}(1,:,:,:)- data{s}(2,:,:,:);%%Extrapolate 
    end                                          %%end if extrapolation
    data{s}(2,:,:,:) = abs(data{s}(2,:,:,:));    %%absolute value of error
    error(1) = max(max(max(data{s}(2,:,:,:))));  %%max step error
    fprintf('Max step error = %e\n',error(1));   %%print max step error
  end                                            %%end if errorchecks
  error(2) =0;                                   %%initial sampling errors
  if r.n.ensemble>1                              %%if ensemble averaging
    data{s}(3,:,:,:)=real(sqrt(data{s}(3,:,:,:)/(r.n.ensemble-1))); %%st.dev.
    error(2) = max(max(max(abs(data{s}(3,:,:,:)))));     %%max of sd matrix
    fprintf('Max sampling error = %e\n',error(2));       %%sampling error
  end                                            %%end if ensemble
  data{s}=reshape(data{s},[3,r.points,r.graphs]);%%reshape space dimensions
end                                              %%end for seqence                                          
if ~strcmp(r.file,'')                            %%if filename not blank
    input = xwrite(error,data,input,raw);        %%save data in file
end                                              %%end if file not blank
if length(input) == 1                            %%if cell array size is 1
    input = input{1};                            %%converts to a structure
end                                              %%end if cell array size is 1
end                                              %%end xsim

%Version 1.03   xsim returns modified input data including default values 


function data1 = xaddcell(data1,data2)
%   data1 = XADDCELL(data1,data2) 
%   Adds two cell arrays together
%   MIT license by Peter D. Drummond, (2015) - see License.txt 
sequence = length(data1);
for s = 1:sequence
    data1{s} = data1{s}+data2{s};
end
end

function raw1 = xrawcell(raw1,raw2)
%   raw1 = XRAWCELL(raw1,raw2) 
%   Combines two cell arrays together
%   MIT license by Peter D. Drummond, (2015) - see License.txt 
sequence = length(data1);
for s = 1:sequence
    data1{s} = data1{s}+data2{s};
end
end
