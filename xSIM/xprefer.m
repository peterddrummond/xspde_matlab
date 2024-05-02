function inlabel = xprefer (in,label,maxl,default,varargin)
%   inlabel  =  XPREFER(in,label,maxl,default,varargin) sets defaults.
%   Input: structure 'in', label 'label', max 'maxl', default `default'.
%   Here 'maxl' is the maximum length of data, set to zero if not needed
%   Output: preferred value for in.label if none is present.
%   The last optional input is for a prior input structure being edited.
%   If input is matrix and default is a cell it is converted to a cell.
%   For partial input of cells or vectors, omitted values use defaults.
%   Truncates input vectors that are greater than 'max' only if 'max'> 0
%   Called by preferences
%   Licensed by Peter D. Drummond, (2024) - see License.txt, XSDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isfield(in,label) && isequal(in.(label),[])   %% If null input present
    inlabel = [];                                %% Leave unchanges
    return;                                      %% Return
end                                              %% End if null input
if ~isfield(in,label)                            %% If no input present
  if nargin > 4 && isfield(varargin{1},label)    %% If any stored input
      inlabel = varargin{1}.(label);             %% Set to previous value
  else                                           %% Else no previous values
      inlabel = default;                         %% Set to default value
  end                                            %% End if prior input
elseif ~iscell(in.(label)) && iscell(default)    %% Default cell, not input
    inlabel = {in.(label)};                      %% Convert to cell
else
    inlabel = in.(label);                        %% use input for label
end                                              %% End if no input present
maxl = max(maxl,length(inlabel));                %% max length --> input
maxl = max(maxl,length(default));                %% max length --> default
if maxl > 0 && ~isequal(default,[])              %% has maxl with default
  while length(default) < maxl                   %% While length too small
    n = max(1,length(default));                  %% get default size
    default(n+1) = default(n);                   %% Set next to default
  end   
  for n = 1:maxl
    if n > length(inlabel) || (iscell(inlabel) && isequal(inlabel(n),{[]}))
      inlabel(n) = default(n);                   %% Set next to default
    end
    if iscell(inlabel) && isvector(inlabel{n})
      while length(inlabel{n}) < length(default{n})
        nv = length(inlabel{n});
        inlabel{n}(nv+1) = default{n}(nv+1);
      end
    end
  end                                            %% End for
end
end                                              %% End function

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XPREFER FUNCTION