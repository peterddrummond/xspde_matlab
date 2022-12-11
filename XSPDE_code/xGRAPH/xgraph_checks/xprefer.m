function inlabel = xprefer (in,label,maxl,default,min,max)
%   inlabel  =  XPREFER(in,label,maxl,default,min,max) sets xGRAPH defaults.
%   Input: structure 'in', label 'label', max length 'maxl', 
%   Default values `default'.
%   Tests out-of-range inputs: minimum value 'min', maximum value 'max'.
%   Here 'maxlen' is the maximum length of data, set to zero if not wanted
%   Output: preferred value for in.label data if none is present.
%   If input is a cell array, it is converted to a matrix.
%   For partial input vectors, fills in any missing values with defaults.
%   Truncates input vectors greater than 'maxl' only if 'maxl'> 0
%   Called by: xpreferences, xgpreferences
%   xGRAPH functions are licensed by Peter D. Drummond, (2021) - see License 
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SET DEFAULTS IF NO INPUT
%
   if ~isfield(in,label)                     %% If no input data present
      inlabel = default;                     %% Set to default value
   else                                      %% Else input data present
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  SET PARTIAL DATA INPUTS
%
      if iscell(in.(label))                  %% If input data is cell
          fprintf ('Warning for input data label "%s"\n',label)
          fprintf ('Converting input cell data to vector type\n');              
          in.(label)=cell2mat(in.(label));   %% Convert to matrix
      end                                    %% End if input data is cell
      sz  = size(in.(label));
      if maxl && maxl <  sz(2)               %% Check max length of data      
        fprintf ('Warning for input data label "%s"\n',label)
        fprintf('Truncating to maximum vector length, size = %d\n', maxl);    
        inlabel(:,1:maxl)=in.(label)(:,1:maxl);%% Truncate in.(label)
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
%  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  TEST INPUT FOR DATA RANGES
%
   sz  = size(inlabel);
   if max>min
       for ind =1:sz(2) 
           if inlabel(1,sz(2))<min
               fprintf ('Warning for input data label "%s"\n',label);
               fprintf ('Minimum value is %d\n',min);
               fprintf('Setting to minimum value, index = %d\n',ind);
               inlabel(1,sz(2))=min;
           end
           if inlabel(1,sz(2))>max
               fprintf ('Warning for input data label "%s"\n',label);
                fprintf ('Maximum value is %d\n',max);
               fprintf('Setting to maximum value, index = %d\n',ind);
               inlabel(1,sz(2))=max;
           end
       end
   end
end                                          %% End function
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END FUNCTION XPREFER
%
