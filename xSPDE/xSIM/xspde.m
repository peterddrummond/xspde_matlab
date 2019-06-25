function [error,input,data,raw] = xspde(input)       
%   [error,input,data,raw] = XSPDE(input) is the xSPDE control function.
%   Input:  input cell array of simulation parameter structures, 'input'.
%   Output: sum of step-size, sampling, comparison errors, labeled 'error',
%   copy of input cell array including preferences, labeled 'input',
%   data cell array of average data, raw data labeled 'data' and 'raw'.
%   licensed by Peter D. Drummond, (2019) - see License.txt 

c = fix(clock);                                  %%Start-up time
fprintf ('\nxSPDE, %d/%d/%d, time %d:%d:%d\n',c);%%Print start-up time
[errorsim,input,data,raw] = xsim(input);         %%Simulation program
errorgraph = xgraph(data,input);                 %%Graphics program   
error = errorgraph+errorsim;                     %%calculate total error
fprintf ('\nxSPDE error %e, time = %f \n\n',error) 
end                                              %%End simulation program