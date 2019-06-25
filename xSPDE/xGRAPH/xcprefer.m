function out = xcprefer (in,label,max,default)
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
if ~isfield(in,label)                        %% If no label data is input
              in.(label) = default;          %% Set in.label to default
end                                          %% End if no label data
if ~iscell(in.(label))                       %% If input data is not cell
         in.(label)={in.(label)};            %% Convert to cell
end                                          %% End if data is not cell
lc = length(in.(label));                     %% Get length of input array
ld = length (default);                       %% Get length of default array
out = cell(1,max);                           %% Initialize output cells
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  CHECK INPUT CELLS, ADD DEFAULT
%                                     
for nc =1:max                                %% Loop over maximum cells
    celldefault = default{min(nc,ld)};       %% Set the default cell
   if nc > lc || isequal(in.(label){nc},[])  %% If count > input or blank
      celldata = celldefault;                %% Set output to  default
   else                                      %% Else input cell exists
      celldata  = in.(label){nc};            %% Set output to input    
      len = length(celldata);                %% Get output length
      for k = 1:length(celldefault)          %% Loop over default length
          if k>len||isequal(celldata(k),{[]})%% If components are missing
          celldata(k)=celldefault(k);        %% Set missing components
          end                                %% End if components missing
      end                                    %% End loop over length
   end                                       %% End if count > input 
   out{nc}=celldata;                         %% Store output cell data
end                                          %% End loop over all cells
end                                          %% End function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END FUNCTION XCPREFER
