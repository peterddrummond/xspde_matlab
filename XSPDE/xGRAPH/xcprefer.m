function inlabel = xcprefer (in,label,max,default)
%   inlabel = XCPREFER (in,label,max,default) sets defaults for cell data.
%   Input: struct 'in', string 'label', max fields 'max', default 'default'.
%   Output:  preferred value for in.label if none is present.
%   If in.(label) is not a cell array, it is converted to a cell array.
%   Here 'default' is a cell array containing vectors of preferences.
%   For missing components or cells, fills in missing values with defaults.
%   If 'default' cell has less fields than 'max', last default cell is used.
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

    if ~isfield(in,label)                    %% If no label data is input
              in.(label) = default;          %% Set to default
    end                                      %% End if no label data
    if ~iscell(in.(label))                   %% If input data is not cell
         in.(label)={in.(label)};            %% Convert to cell
    end                                      %% End if data is not cell
    lc = length (in.(label));                %% Get length of cell array
    for nc = 1:lc                            %% For cells present
        cdata = in.(label){nc};              %% Get input cell data at nc
        nd = min(nc,length(default));        %% Cell index for default data
        while ~ischar(cdata) && (length(cdata) < length(default{nd}))
          ns=length(cdata);                  %% Get length
          if iscell(default{nd})             %% If data is nested cell
            cdata{ns+1} = default{nd}{ns+1}; %% Add cell default at ns+1
          else                               %% If data in cell is not cell
            cdata(ns+1) = default{nd}(ns+1); %% Add matrix default at ns+1
          end                                %% End if is nested cell
        end                                  %% End while
        inlabel{nc} = cdata;                 %% Store current cell
    end                                      %% End for cells present
    for nc =lc+1:max                         %% For missing cells   
        inlabel{nc} = default{min(nc,length(default))};  %% Set default
    end                                      %% End for missing cells
end                                          %% End function