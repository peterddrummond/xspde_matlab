function [error,input,data,raw] = xsim(input)                 
%   [error,input,data,raw] = XSIM(input) 
%   Solves stochastic differential and partial differential equations.
%   It returns averages as defined by the 'input' cell array. 
%   XSIM also calls parallel loops, estimates errors, writes files.
%   The 'error' is the sum of max sampling and step errors over sequences. 
%   The output cell array 'input' is input data including default values.
%   The output cell array 'data' includes the averages and error bars.
%   The output cell array 'raw' is the raw trajectories if required.
%   Licensed by Peter D. Drummond, (2015) - see License.txt 

tic();                                           %%set timer

input = xmakecell(input);                        %%make input a cell
rawdata =    xprefer(input{1},'rawdata',0,'');
rawinput = {};                                   %%initialize raw input
if rawdata ~=  ''                                %%If data is character
    [rawinput,~,raw] = xread(rawdata);           %%Input rawdata from file 
end                                              %%End if rawdata 
rawinput = xmakecell(rawinput);                  %%change rawinput to cell
input = xpreferences(input,rawinput);            %%get 'input' defaults
sequence = length(input);                        %%get sequence length
data = cell(sequence);                           %%allocate data as cell
l = cell(sequence);                              %%allocate lattice as cell

%%Loop over sequence to initialize data sets

for s=1:sequence                                 %%loop over sequence
  l{s} = xlattice(input{s});                     %%initialise lattice
  if input{s}.seed >=0  && ~input{s}.octave      %%check if not octave
    sd = RandStream('CombRecursive','Seed',1);   %%define RNG type
    RandStream.setGlobalStream(sd);              %%set RNG type
  end
  for n=1:l{s}.functions
        data{s}{n} = zeros(l{s}.d.data{n});  %%initialise data cells
  end
  fprintf ('%s sequence %d: %s\n',l{s}.version,s,l{s}.name);%%version name
  if l{s}.print >1                               %%if print switch verbose
    fprintf ('\n xSIM parameters\n');            %%display input data
    display (l{s});                              %%display input data
  end                                            %%end if print switch
end                                              %%end loop over sequence
r = l{1};                                        %%first sequence cell 

%%Carry out parallel simulations if needed

if r.ensembles(3) > 1                            %%if parallel simulations
    parfor npar = 1:r.ensembles(3)               %%loop on parallel threads
        [dp,raw1] = xensemble(npar,l);           %%call xensemble function
        raw(:,:,:,npar) =raw1;                   %%accumulate raw data
        data = xaddcell(data,dp);                %%accumulate averages
    end                                          %%end parallels loop
    raw = reshape(raw,[sequence,r.errorchecks,r.ncopies]); %%reshape raw 
else                                             %%no parallels specified
    [data,raw] = xensemble(1,l);                 %%call xensemble function
end                                              %%end if parallels

%%Loop over sequence to finalize data sets

error = 0;                                       %%initialize sum of errors
for s=1:sequence                                 %%loop over sequence
  r = l{s};                                      %%get lattice for sequence
  esp =0;                                        %%initial sampling errors
  es =0;                                         %%initial step errors
  for n=1:r.functions
    data{s}{n}(:,3,:)=data{s}{n}(:,3,:) - abs(data{s}{n}(:,1,:)).^2;
    if r.ncopies > 1                               %%if ensemble averaging   
      data{s}{n}(:,3,:)=real(sqrt(data{s}{n}(:,3,:)/(r.ncopies-1))); 
      esp = esp+max(max(max(abs(data{s}{n}(:,3,:))))); %%max of sd matrix
    end                                          %%end if ensemble
    if r.errorchecks > 1                         %%if errorchecks needed
      data{s}{n}(:,2,:)=data{s}{n}(:,2,:) - data{s}{n}(:,1,:);
      if r.order > 0                             %%if extrapolation order
        data{s}{n}(:,2,:) = data{s}{n}(:,2,:)/(2^(r.order)-1.);
        data{s}{n}(:,1,:) = data{s}{n}(:,1,:)- data{s}{n}(:,2,:);
      end                                        %%end if extrapolation
      data{s}{n}(:,2,:) = abs(data{s}{n}(:,2,:));
      es = es+ max(max(max(data{s}{n}(:,2,:)))); %%Sum of max step errors
    end                                          %%end if errorchecks
    data{s}{n}=reshape(data{s}{n},r.gpoints{n});
  end
  if esp>0      
      fprintf('Sum of max sampling errors, sequence %d = %e\n',s,esp);
  end
  if es>0    
      fprintf('Sum of max step errors, sequence %d = %e\n',s,es);
  end
  error = error+esp+es;                          %%store sum of errors
end                                              %%end seqence  

%%Write data sets to output file

if ~strcmp(r.file,'')                            %%if filename not blank
    input = xwrite(error,data,input,raw);        %%save data in file
end                                              %%end if file not blank
if length(input) == 1                            %%if  input length is 1
    input = input{1};                            %%returns struct, not cell
end                                              %%end if input length is 1
timer = toc();
fprintf('xSIM sequence completed, time = %f \n\n',timer);%%print time taken
error = [error,timer];
end                                              %%end xsim



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
