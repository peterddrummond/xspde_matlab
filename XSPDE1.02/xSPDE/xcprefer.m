function inlabel = xcprefer (in,label,max,default)
%   inlabel = XCPREFER (in,label,max,default) sets defaults for cell data.
%   Input: struct 'in', string 'label', max fields 'max', default 'default'.
%   Output:  preferred value for in.label if none is present.
%   If in.(label) is not a cell array, it is converted to a cell array.
%   Here 'default' is a cell array containing vectors of preferences.
%   For missing components or cells, fills in missing values with defaults.
%   If 'default' cell has less fields than 'max', last default cell is used.
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

    if ~isfield(in,label)                    %% If no label data
              in.(label) = default;          %% Set to default
    end                                      %% End if no label data
    if ~iscell(in.(label))                   %% If label data not cell
         in.(label)={in.(label)};            %% Convert to cell
    end                                      %% End if label data not cell
    lc = length (in.(label));                %% Get length of cell array
    for nc =1:lc                             %% For cells present
        nd = min(nc,length(default));        %% Index for default cell data
        while length(in.(label){nc}) < length(default{nd}) %% While short
          ns=length(in.(label){nc});                       %% Get length
          in.(label){nc}(ns+1) = default{nd}(ns+1);        %% Add default
        end                                                %% End while
    end                                      %% End for cells present
    for nc =lc+1:max                         %% For missing cells   
       in.(label){nc} = default{min(nc,length(default))};  %% Set default
    end                                      %% End for missing cells
    inlabel=in.(label);                      %% Return modified inlabel
end                                          %% End function