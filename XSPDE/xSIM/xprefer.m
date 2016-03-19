function inlabel = xprefer (in,label,max,default)
%   inlabel  =  XPREFER(in,label,max,default) sets defaults for input data.
%   Input: structure 'in', label 'label', max 'max', default `default'.
%   Output: preferred value for in.label if none is present.
%   If input is a cell array, it is converted to a matrix.
%   For partial input vectors, fills in any missing values with defaults.
%   Truncates input vectors that are greater than 'max' only if 'max'> 0
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 

   if ~isfield(in,label)                    %% If no label present
      in.(label) = default;                  %% Set default value
    else                                     %% Else label present
      if iscell(in.(label))                  %% If label data is cell
         in.(label)=cell2mat(in.(label));    %% Convert to matrix
      end                                    %% End if label data is cell
      while length(in.(label)) < length(default)  %% While length too small
        n=length(in.(label));                %% get current size
        in.(label)(n+1) = default(n+1);      %% Set next to default
      end                                    %% End while
    end                                      %% End if no label
    if max && max < length(in.(label))
      inlabel(1:max)=in.(label)(1:max);      %% Truncate in.(label)
    else
      inlabel=in.(label);                    %% Return modified in.(label)
    end
end                                          %% End function
