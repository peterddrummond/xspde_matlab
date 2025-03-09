function e1 = GBS100( )
%Batch testing script for xqsim program
%Matrix size      = 100*100
%Sequence length  = 12
%Total tests      = 68
%Total graphs     = 152
%Chi-square error = 1.0 +/- 0.2
%Typical timing   = 700s

p.matrix     = @Identity;                        %matrix type
p.dimensions = 0;                                %dimension
p.phase      = 1;
p.modes      = 100;                              %matrix size m
p.verbose    = 1;
p.name       = sprintf('+P non-uniform GBS, M=%d',p.modes); 
p.tr         = 0.5*ones(1,p.modes);              %transmission
I            = ones(1,p.modes/5);                %identity vector
p.sqz        = [I/4,2*I,4*I,I/2,I];              %nonuniform squeezing
p.thermal    = zeros(1,p.modes);                 %decoherence factor 
p.correl{4}  = 1:5;                              %correlation order 
p.correl{6}  = 1:5;                              %click order 
p.part{7}    = p.modes;                          %One-fold partititon 
p.part{8}    = p.modes/2*[1,1];                  %Two-fold partititon    
p.part{9}    = p.modes/4*[1,1,1,1];              %Four-fold partition
p.ensembles  = [100,10,12];                      %repeats for errors
p.cutoff     = 1.e-7;
p.observe    = {@Np,@X2,@Y2,@Nm,@K,@Km,@K1,@Kn,@Kn};
p.compare    = {@Npc,@X2c,@Y2c,@Nmc,@Kc,@Kmc,@K1c,@Knc,@Knc};
p.glabels    = {{{},'Mode j'},{{},'Mode j'},{{},'Mode j'},{{},'Order'},...
               {{},'Mode j'},{{},'Clicks m'},{{},'Clicks m_1','Clicks m_2'},...
               {{},'Clicks m_1','Clicks m_2','Clicks m_3','Clicks m_4'}};
p.olabels    = {'<n>','<x^2>','<y^2>','<n_1...n_j>','<G^{(1)}>',...
               '<G^{(n)}>','G_1(m)','G_2(m)','G_4(m)'};
p.diffplot   = {2,2,2,2,2,2,2,2,2};
p.xk{4}      = p.correl(4);
p.xk{6}      = p.correl(6);
p.xk{7}      = {0:p.modes};
p.xk{8}      = {0:p.modes/2,0:p.modes/2};
p.xk{9}      = {0:p.modes/4,0:p.modes/4,0:p.modes/4,0:p.modes/4};
p0 = p;

p.transfer   = @Initialgaussian;
p1           = p;
p1.thermal   = 0.5*ones(1,p.modes);                  %decoherence factor
p1.name      = sprintf('+P GBS thermalized, M=%d',p.modes); %test name
p2           = p;
p2.sqz       = ones(1,p.modes);
p2.thermal   = zeros(1,p.modes);                     %decoherence factor
p2.name      = sprintf('+P uniform GBS, M=%d',p.modes);     %test name

p3           = p;
p3.sqz       = ones(1,p.modes);
p3.thermal   = ones(1,p.modes);                      %decoherence factor
p3.matrix    = @Unitary;                             %random unitary
p3.name      = sprintf('+P  unitary thermal, M=%d',p.modes);  %test name

pW           = p;
pW.averages  = 1:4;
pW.phase     = 2; 
pW.name      = sprintf('W non-uniform GBS, M=%d',p.modes);  %test name

p1W          = p1;
p1W.averages = 1:4;
p1W.phase    = 2;                              
p1W.name     = sprintf('W non-uniform thermalized, M=%d',p.modes);%test name

p2W          = p2;
p2W.averages = 1:4;
p2W.phase    = 2;                                 
p2W.name     = sprintf('W uniform GBS, M=%d',p.modes);      %test name

p3W          = p3;
p3W.averages = 1:4;
p3W.phase    = 2;                              
p3W.name     = sprintf('W unitary thermal, M=%d',p.modes);       %test name

pQ           = pW;
pQ.phase     = 3; 
pQ.name      = sprintf('Q non-uniform GBS, M=%d',p.modes);  %test name

p1Q          = p1W;
p1Q.phase    = 3;                              
p1Q.name     = sprintf('Q non-uniform thermalized, M=%d',p.modes);%test name

p2Q          = p2W;
p2Q.phase    = 3;                                 
p2Q.name     = sprintf('Q uniform squeezed, M=%d',p.modes);      %test name

p3Q          = p3W;
p3Q.phase    = 3;                              
p3Q.name     = sprintf('Q unitary thermal, M=%d',p.modes);       %test name

e1    = xspde(p0,p1,p2,p3,pW,p1W,p2W,p3W,pQ,p1Q,p2Q,p3Q);