function [a,av,avp,raw] = xpathfb(v,nc,p) 
%   [a,av,avp,raw] = XPATHFB(v,nc,p)  solves FBSDEs by Picard iteration.
%   v is a cell of initial randoms, nc the check index, p are parameters
%   Fields are a ={x,y}; x propagates forward in time, y backward
%        nc is the check index: nc = 1 for fine, nc = 2 for coarse
%        If nc = 1, and p.checks = 1, the checks loop is taken twice
%   Input:  initial random values v, check index nc, parameters p
%   Output: Final field cells a, averages av, RMS iteration error e
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Called by: xensemble
%   Calls: firstfb, initial, deriv, xnoisefb, xdatafb
%   Needs: d.raw, d.av, fieldcells, fbcells, origins, ranges, dtr,
%          points, steps, checks, fields, ensembles, propagator
%   Licensed by Peter D. Drummond, (2025) - see License.txt, XSPDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

raw = czeros(p.d.raw);                           % initialize raw cells
p.f = 1:p.fieldcells;
pts = p.points{1}(1);  
a = cell(1,p.fbcells);
a1 = a; a2 = a; b = a;
dt = p.dtr;
sgn   = sign(p.iterfb);
iters = abs(p.iterfb);
p.b = (1+p.fieldcells):p.fbcells;
p.fbsteps = p.steps(1)*(1+(2-nc)*p.checks);     % increase the local steps                                            
p.fbpoints = 1+(p.points{1}-1)*p.fbsteps;
wstore = xnoisefb(p);                            % compute stored noise
for c = 1:p.fbcells  
  p.fbsize{c} = [p.fields{c},p.fbpoints,p.ensembles(1)];
  a1{c} = p.firstfb{c}(p)+zeros(p.fbsize{c});    % set a to first value
end 
p.t = p.origins(1)+p.dtr*(0:(p.fbpoints-1));
av  = czeros(p.d.av);                            % initialize average
for c = 1:p.fbcells 
  a1{c} = a1{c}+zeros(p.fbsize{c});
end 
for iter = 1:iters                               % loop over iterations      
  for c = 1:p.fbcells 
    a{c} = reshape(a1{c}(:,1,:),p.d.ca{c});
  end   
  n = 1+(p.fbpoints-1)*(sgn<0);                  % Initialize time index
  npnoise = 0;                                   % Initialize noise index
  p.t = p.origins(1)+(p.ranges(1))*(sgn<0);      % Initial time
  for fbind = (1:-2:-1)*sgn                      % Loop on forward/backward
    ci = (fbind>0)*p.f + (fbind<0)*p.b;          % Directional indices
    cr = (fbind>0)*p.b + (fbind<0)*p.f;          % Reversed indices
    inpt = 1+(fbind<0)*(p.fbpoints-1);           % Initial point index 
    for np = 1:pts                               % loop until time tmax
      if np == 1                                 % If first point
        for c = ci                               % Loop on current cells
          a{c} = zeros(p.d.ca{c});               % Set cells to zeros
          a{c} = a{c}+p.initial{c}(a{cr},v{:},p);% Initialise x,y fields
          a2{c}(:,inpt,:) = a{c};                % Initialize total field
          npnoise = npnoise+fbind*sgn;           % increment noise index
        end                                      % end loop on current
      else                                       % else not first point
        for step = 1:p.steps(1)                  % loop over steps
          for check = 1:(1+p.checks(1))          % loop over checks          
            z = wstore{npnoise};                 % get  noise;     
            npnoise = npnoise+fbind*sgn;         % increment noise index
            if nc == 2                           % if low res,
              if check == 1                      % if low res, first check
                zstore = z;                      % store current noise
                continue                         % exit the check loop
              else                               % else low res, 2nd check
                for c=1:p.noisecells(1)          % loop over all noises
                  z{c} = 0.5*(z{c}+zstore{c});   % average noise
                end
              end                                % end if first check 
            end                                  % End if low res     
            n1 = n+fbind;
            for c = 1:p.fbcells                  % loop over all cells
              b{c} = reshape((a1{c}(:,n,:)+a1{c}(:,n1,:))/2,p.d.ca{c});  
            end
            for c = ci
              a{c} = p.propagator{c}.*a{c};      % propagate 
              a{c} = a{c} + reshape(p.deriv{c}(b{:},z{:},p)*dt,p.d.ca{c});
              a{c} = p.propagator{c}.*a{c};      % propagate 
              a2{c}(:,n1,:) = a{c}; 
            end
            n = n1;
            p.t = p.t+p.dtr*fbind;               % Increment time
         end                                     % End check loop
        end                                      % End for steps
      end                                        % End if first point
    end                                          % End time loop
  end                                            % End fb loop
  a1 = a2;
  avp = av;
  [av,delt] = xdatafb(a2,av,p);                  % Compute all averages
  if p.verbose > 1  && iters > 1                 % if verbose printing
   fprintf('FB convergence error = %d, iteration = %d\n',delt,iter);   
   fprintf('Typical end values = {%d,%d}',a2{1}(1,pts,1),a2{2}(1,1,1));
  end                                            % end if verbose printing   
end                                              % End iteration loop
if p.rawdata
    raw = a2;
end
if p.verbose > 0  && iters > 1                   % if verbose printing
   fprintf('FB convergence error = %d, iteration = %d\n',delt,iter);
end                                              % end if verbose printing   
end                                              % end trajectory function