function [error,data,input,raw] = xspde(input)       
%   [error,input,data,raw] = XSPDE(input) is the xSPDE control function.
%   Input:  input cell array of simulation parameter structures, 'input'.
%   Output: 'error' vector, 'data' cell array of average data,           
%          'input' cell array with parameters,'raw' data of trajectories.         
%   licensed by Peter D. Drummond, (2022) - see License.txt
c = fix(clock);                                    %%Start-up time
fprintf('\n\nxSPDE, %d/%d/%d, start-up time %d:%d:%d\n',c);%%Start-up time
[error,data,input,raw] = xsim(input);              %%Simulation program
xgraph(data,input);                                %%Graphics program 
end                                                %%End simulation program