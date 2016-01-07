function [error,input,data] = xspde(input)       
%   [error,input,data,] = XSPDE(input,graph) is the xSPDE control function.
%   Input:  input cell array, 'input'.
%   Output: summed errors 'error' of step-size, sampling, comparison errors
%   Output: input cell array including preferences, labeled 'input'
%   Output: data cell array of average data, labeled 'data'
%   MIT licensed by Peter D. Drummond, (2015) - see License.txt 

fprintf ('\nxSPDE starting\n');                  %%integration number
[errorsim,input,data] = xsim(input);             %%Simulation program
errorgraph = xgraph(input,data);                 %%Graphics program   
error = errorgraph+errorsim;                     %%calculate max error

end                                              %%End simulation program

%Version 1.04   Simulation functions return preferred input defaults