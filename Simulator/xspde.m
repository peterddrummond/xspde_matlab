function [error,data,input,raw] = xspde(varargin)       
%   [error,input,data,raw] = XSPDE(input) is the xSPDE control function.
%   Input:  list of simulation parameter structures.
%   Output: 'error', 'data' average, 'input' parameters, 'raw' trajectories.         
%   Licensed by Peter D. Drummond, (2023) - see License.txt, XSPDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[error,data,input,raw] = xsim(varargin{:});      %%Simulation program
xgraph(data,input);                              %%Graphics program 
end                                              %%End XSDE program

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XSDE FUNCTION