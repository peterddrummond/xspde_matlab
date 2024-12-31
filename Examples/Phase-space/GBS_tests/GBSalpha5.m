function e1 = GBSalpha5( )

p.dimensions = 0;
p.phase      = 1;                                %+P phase-space
p.verbose = 1;
p.modes      = 5;                                %matrix size m
p.name       = sprintf('+P coherent, M=%d',p.modes); 
p.tr         = .5*ones(1,p.modes);               %transmission
I            = ones(1,p.modes/5);                %identity vector
p.sqz        = [I/2,I,1.5*I,2*I,0*I];            %nonuniform squeezing
p.alpha      = [I/4,2*I,4*I,I/2,I];              %nonuniform squeezing
p.thermal    = 0.5*ones(1,p.modes);              %decoherence factor 
p.ensembles  = [1000,10,1];                       %repeats for errors
p.observe    = @pn;
p.compare    = @nc;
p.glabels    = {{' ','Mode j'}};
p.olabels    = {'<n>'};
p.diffplot   = {1};
e1           = xspde(p);