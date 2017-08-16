function inlabel = xprefer (in,label,max,default)
%   inlabel  =  XPREFER(in,label,max,default) sets defaults for vector data.
%   Input: structure 'in', label 'label', max 'max', default `default'.
%   Here 'max' is the maximum length of data, set to zero if not wanted
%   Output: preferred value for in.label if none is present.
%   If input is a cell array, it is converted to a matrix.
%   For partial input vectors, fills in any missing values with defaults.
%   Truncates input vectors that are greater than 'max' only if 'max'> 0
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License 

   if ~isfield(in,label)                     %% If no input data present
      inlabel = default;                     %% Set to default value
   else                                      %% Else input data present
      if iscell(in.(label))                  %% If input data is cell
          fprintf ('Warning for input data label "%s"\n',label)
          fprintf ('Converting input cell data to vector type\n');              
          in.(label)=cell2mat(in.(label));   %% Convert to matrix
      end                                    %% End if input data is cell
      sz  = size(in.(label));
      if max && max <  sz(2)                 %% Check max length of data      
        fprintf ('Warning for input data label "%s"\n',label)
        fprintf('Truncating to maximum vector length, size = %d\n', max);    
        inlabel(:,1:max)=in.(label)(:,1:max);%% Truncate in.(label)
        sz  = size(inlabel);
      else                                   %% Else length of data is OK
        inlabel=in.(label);                  %% Return modified in.(label)
      end                                    %% End check length of data
      szd = size(default);
      while sz(2) < szd(2)                   %% While length too small
        sz(2) =   sz(2)+1;
        inlabel(:,sz(2)) = default(:,sz(2));%% Set next to default
      end                                    %% End while   
   end                                       %% End if no input data 
end                                          %% End function
