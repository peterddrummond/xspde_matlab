function p = phasepreferences(p)
% p = phasepreferences(p) sets phase-space xpreferences for xspde
% Input and output are parameter structures
% p.phase = 1,2,3 for normal, symmetric, and anti-normal ordering 
% Licensed by Peter D. Drummond, (2023) - see License.txt 

% INITIALIZE STANDARD INPUTS

p.modes   = xprefer(p,'modes',1,1);              %% default mode number
p.samples = xprefer(p,'samples',1,1);            %% default sample number
p.minvar  = xprefer(p,'minvar',1,10^-100);       %% set minimum variance
p.thermal = xprefer(p,'thermal',1,0);            %% Default thermal 
p.tr      = xprefer(p,'tr',1,1);                 %% Default transmission
p.pnames  = xprefer(p,'pnames',0,{'+P','W','Q'});%% phase-space names
p.gauge   = xprefer(p,'gauge',1,0);              %% phase-space gauge
p.counts  = xprefer(p,'counts',1,0);             %% simulates counts 
p.jump    = xprefer(p,'jump',1,0);               %% phase-space jump switch
if p.counts == 1
  p.path     = xprefer(p,'path',0,@countpath);
  p.jmethod  = xprefer(p,'jmethod',0,@phasejump);
  p.uniforms = p.modes;                          %% required for counts
  p.gauge    = 1;                                %% required for counts
end
p.initial =   xprefer(p,'initial',0,{@Initialgaussian});
if (isfield(p,'matrix'))
  p.Tname   = xprefer(p,'Tname',0,p.matrix(0));  %%Default matrix name
else
  p.matrix  = xprefer(p,'matrix',0,@(p) eye(p.modes)); %%Default matrix 
  p.Tname   = 'Identity';
end
p.inrandoms = xprefer(p,'inrandoms',1,2*p.modes*(min(p.phase,2)));
p.sqz     =   xprefer(p,'sqz',1,0);                %%Default squeezing
p.sqz     =   [p.sqz,0*(1:(p.modes-length(p.sqz)))];
p.alpha   =   xprefer(p,'alpha',1,0);              %%Default amplitude
p.alpha   =   [p.alpha,0*(1:(p.modes-length(p.alpha)))];
p.correl  =   xprefer(p,'correl',0,{[]});          %% correlations 
p.part    =   xprefer(p,'part',0,{[]});            %% partitions 
p.pname   =   p.pnames{p.phase};                   %% current name
if p.phase == 1 
    p.fields = xprefer(p,'fields',1,{p.gauge+2*p.modes(1)});
else
    p.fields = xprefer(p,'fields',1,{p.modes(1)});
end
if p.jump == 1
     p.jmethod =  xprefer(p,'jmethod',0,@phasejump);
     p.fields{2} = p.modes;
     p.deriv = xprefer(p,'deriv',1,{@(a,~,~,~,~) 0*a});
     p.deriv{2} = @(~,b,~,~,~) 0*b;
     p.initial{2} = @(~,~,~) zeros([p.modes,1]);
end
end                                              %%end function