function out = xcprefer (in,label,maxc,default)
%   out = XCPREFER (in,label,max,default) sets defaults for cell data.
%   Input: struct 'in', string 'label', max fields 'max', default 'default'
%   Output:  preferred value for in.label if none is present.
%   If in.(label) is not a cell array, it is converted to a cell array.
%   Here 'default' is a cell array containing cells/vectors of preferences.
%   For missing components or cells, fills in missing values with defaults.
%   If 'default' has less fields than 'max', last default cell is used.
%   xSPDE functions are licensed by Peter D. Drummond, (2019) - see License 
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CHECK INPUT DATA TYPES
%
cells = iscell(default);                     %% Set flag for a default cell                                   
if cells                                     %% If default is cell type
    maxv = length(default{1});               %% Set component numbers
else                                         %% If default not cell type
    maxv =  maxc;                             %% Set components to max
    default = {default};                     %% Set default to cell
    maxc =1;                                  %% Set max cells to 1
end                                          %% End if cell data
if ~isfield(in,label)                        %% If no label data is input
              in.(label) = default;          %% Set in.label to default
end                                          %% End if no label data
if ~iscell(in.(label))                       %% If input data is not cell
         in.(label)={in.(label)};            %% Convert to cell
end                                          %% End if data is not cell
lc = length(in.(label));                     %% Get length of input array
ld = length (default);                       %% Get length of default array
out = cell(1,maxc);                           %% Initialize output cells
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CHECK INPUT CELLS, ADD DEFAULT
%                                     
for nc =1:maxc                                %% Loop over maximum cells
    cdefault = default{min(nc,ld)};          %% Set the default cell
   if nc > lc || isequal(in.(label){nc},[])  %% If count > input or blank
      cdata = cdefault;                      %% Set output to  default
   else                                      %% Else input cell exists
      cdata  = in.(label){nc};               %% Set output to input    
      len = length(cdata);                   %% Get input vector length
      dlen = length(cdefault);               %% Get default vector length
      for k = 1:max(dlen,maxv)               %% Loop over output length
          if k>len||isequal(cdata(k),{[]})   %% If components are missing
          cdata(k)=cdefault(min(k,dlen));    %% Set the default component
          end                                %% End if components missing
      end                                    %% End loop over length
   end                                       %% End if count > input 
   if cells                                  %% If default is cell type
       out{nc}=cdata;                        %% Store output cell data
   else                                      %% Else default isnt cell type
       out=cdata;                            %% Store output non-cell data
   end                                       %% End if default is cell type
end                                          %% End loop over all cells
end                                          %% End function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END FUNCTION XCPREFER
