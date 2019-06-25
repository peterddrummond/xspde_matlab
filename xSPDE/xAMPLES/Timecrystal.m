function [e] = Timecrystal()                            
%   e  =  Timecrystal() tests xSPDE for a subharmonic quantum oscillator
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License                         

in.name =       'Timecrystal';
in.ensembles =  [1,1];                       %%samples,ensembles,parallel
in.initial   =  @(w,r)    0+0*w ;                %%Initialisation  handle
in.steps     =  200;
in.step     =  @xRK4;
in.ranges   =   100;
in.fields   =   2;    
in.C         =  1;
in.L         =  4;
in.da        =  @Timeda;                        %%Derivative  handle
in.observe{1} = @(a,r) (a(1,:)+a(2,:))/2;
in.observe{2} = @(a,r) a(1,:).*a(2,:);
in.olabels    = {'a','n'};
in.xlabels =    {'\tau'};
e         =     xspde(in);                       %%Stochasic program
end

function da=Timeda(a,xi,r)
a2=(1-a.*a);
a2s=sqrt(a2);
b=flip(a,1);
da=-(r.C+0.5)*a+r.L*b.*a2+a2s.*xi;
end