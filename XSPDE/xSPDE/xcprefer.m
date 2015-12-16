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
        cdata = in.(label){nc};              %% Get cell data
        while iscell(cdata)                  %% While data is nested cell 
            cdata = cdata{1};                %% Un-nest nested cell 
        end                                  %% End if data is nested cell 
        nd = min(nc,length(default));        %% Cell index for default data
        while length(cdata) < length(default{nd}) %% While too short
          ns=length(cdata);                  %% Get length
          cdata(ns+1) = default{nd}(ns+1);   %% Add default
        end                                  %% End while
        in.(label){nc} = cdata;              %% Store cell
    end                                      %% End for cells present
    for nc =lc+1:max                         %% For missing cells   
       in.(label){nc} = default{min(nc,length(default))};  %% Set default
    end                                      %% End for missing cells
    inlabel=in.(label);                      %% Return modified inlabel
end                                          %% End function

%   Version 1.03 - removes nested cells if present in the input data