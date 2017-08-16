function osize = xfcheck (fname,fnumber,oarray,esize)
%   inlabel  =  XFCHECK(fname,fnumber,oarray,esize) checks function data.
%   Input: function name 'fname', 
%          function number 'fnumber', 
%          output 'oarray', 
%          expected array size 'esize'.
%   Throws error if user function returns the wrong array size
%   If first expected size is zero, only checks the other array sizes
%   xfcheck returns the actual array size found
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

 osize = size(oarray);

 if fnumber > 0
     fn = sprintf('{%d}',fnumber);
 else
     fn = '';
 end
 if esize(1)>0
   if ~isequal(osize,esize)
     fprintf('User defined "%s%s" returned a %d ',fname,fn,osize(1));
     fprintf('x %d ', osize(2:end)); fprintf('array\n');
     fprintf('Should return a %d ',esize(1));
     fprintf('x %d ', esize(2:end)); fprintf('array\n');
     error('%s%s returns wrong array size',fname,fn);
   end
 else
   if ~isequal(osize(2:end),esize(2:end))
     fprintf('User defined "%s%s" returned a %d ',fname,fn,osize(1));
     fprintf('x %d ', osize(2:end)); fprintf('array\n');
     fprintf('Should return a (*) ');
     fprintf('x %d ', esize(2:end)); fprintf('array\n');
     error('%s%s returns wrong array size',fname,fn);
   end
end

