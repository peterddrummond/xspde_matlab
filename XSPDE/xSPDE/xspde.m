function [error,input,data] = xspde(input)       
%   [error,input,data,] = XSPDE(input,graph) is the xSPDE control function.
%   Input:  input cell array, 'input'.
%   Output: error vector 'error' with step-size, sampling, comparison errors
%   Output: input cell array including preferences, labeled 'input'
%   Output: data cell array of average data, labeled 'data'
%   MIT licensed by Peter D. Drummond, (2015) - see License.txt 

[error,input,data] = xsim(input);                %%Simulation program
error(3) = xgraph(input,data);                   %%Graphics program   
error = max(error);                              %%calculate maximum error
end                                              %%End simulation program

%Version 1.03   Simulation functions return preferred input defaults