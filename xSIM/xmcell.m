function input = xmcell (input)
%   incell = XMCELL (input) makes cells from non-cell data.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

    if ~iscell(input) &&  ~isequal(input,[]) %% If input not cell or void
         input = {input};                    %% Convert to cell
    end                                      %% End if input not cell
end                                          %% End function