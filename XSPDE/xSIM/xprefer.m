function inlabel = xprefer (in,label,max,default)
%   inlabel  =  XPREFER(in,label,max,default) sets defaults for input data.
%   Input: structure 'in', label 'label', max 'max', default `default'.
%   Here 'max' is the maximum length of data, set to zero if not wanted
%   Output: preferred value for in.label if none is present.
%   If input is a cell array, it is converted to a matrix.
%   For partial input vectors, fills in any missing values with defaults.
%   Truncates input vectors that are greater than 'max' only if 'max'> 0
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 

   if ~isfield(in,label)                    %% If no input data present
      in.(label) = default;                  %% Set to default value
    else                                     %% Else input data present
      if iscell(in.(label))                  %% If input data is cell
         in.(label)=cell2mat(in.(label));    %% Convert to matrix
      end                                    %% End if input data is cell
      while length(in.(label)) < length(default)  %% While length too small
        n=length(in.(label));                %% get current size
        in.(label)(n+1) = default(n+1);      %% Set next to default
      end                                    %% End while
    end                                      %% End if no input data 
    if max && max < length(in.(label))       %% Check max length of data 
      inlabel(1:max)=in.(label)(1:max);      %% Truncate in.(label)
    else                                     %% Else length of data is OK
      inlabel=in.(label);                    %% Return modified in.(label)
    end                                      %% End check length of data
end                                          %% End function
