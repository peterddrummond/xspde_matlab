function [] = xpr(m,p,varargin)
%   inlabel  =  XPRINT(sw,p,varargin) prints xSPDE messages
%   Prints output if m <= p.print.
%   Minimal if p.print = -1       Prints start-up and faults
%   Brief if   p.print =  0:     +Prints total integration errors
%   Normal if  p.print =  1:     +Prints function errors and progress
%   Verbose if p.print =  2:     +Prints internal parameters as well
%   
%   xSPDE functions licensed by Peter D. Drummond, (2021) - see License  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  if m <= p.print
    fprintf (varargin{:});
  end
end      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END FUNCTION XPR
%
