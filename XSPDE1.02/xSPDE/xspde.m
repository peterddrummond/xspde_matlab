function [e,data] = xspde(input)       
%   [e,data] = XSPDE(input,graph) is the main xSPDE control function.
%   Input:  input cell array, 'input'.
%   Output: error vector 'e' with step-size, sampling, comparison errors
%   Output: data cell array of average data, labeled 'data'
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

[e,data] = xsim(input);                      %%Simulation program
e(3) = xgraph(data,input);                   %%Graphics program   
e = max(e);                                  %%calculate maximum error
end                                          %%End simulation program
