function incell = xmakecell (input)
%   incell = XMAKECELL (input) makes cells from non-cell data.
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

    if ~iscell(input)                        %% If input data not cell
         incell ={input};                    %% Convert to cell
    else                                     %% Else already cell
        incell =input;                       %% Don't convert to cell
    end                                      %% End if input not cel
end                                          %% End function