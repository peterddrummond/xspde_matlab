
function inlabel = xprefer (in,label,default)
%   inlabel  =  XPREFER(in,label,default) sets defaults for matrix data.
%   Input:   structure 'in', label 'label',  default value `default'.
%   Output:  preferred value for in.label if none is present.
%   If input is a cell array, it is converted to a matrix.
%   For partial input vectors, fills in any missing values with defaults.
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

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
    inlabel=in.(label);                      %% Return modified in.(label)
end                                          %% End function
