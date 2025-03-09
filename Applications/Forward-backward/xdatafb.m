function [av,e] = xdatafb(ac,avp,p) 
%   [av,del] = XDATAFB(ac,avp,p)  stores data from a 
%   simulation of a  forward-backward stochastic differential equation
%   Used for solving FBSDEs & FBPSDEs by iteration.
%   Calls xdata for each time point and required output
%   Input:  trajectory cells ac, previous averages avp, parameters p
%   Output: new cell array of averages av, RMS iteration error e
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Called by: xpathfb
%   Calls:     xdata
%   Needs: d.av, fbcells, origins, dt, fbcells, d.ca, 
%          averages, noisefactor, fbpoints, fbsteps
%   Licensed by Peter D. Drummond, (2025) - see License.txt, XSPDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%                       

npt = 0;                                         % initialise index count
e = 0;                                           % initialise RMS error
av = czeros(p.d.av);                             % initialise averages
a  = cell(1,p.fbcells);                          % initialise local fields
for np = 1:p.fbsteps:p.fbpoints                  % loop over all indices
  npt = 1+npt;                                   % update index count
  p.t = p.origins(1)+(npt-1)*p.dt;               % compute time 
  for c = 1:p.fbcells                            % loop over all cells
    a{c} = reshape(ac{c}(:,np,:),p.d.ca{c});     % current field forward
  end                                            % end loop over all cells
  for n = p.averages                             % loop over averages  
      av{n}(:,npt,:) = xdata(a,n,p);             % store time-domain data
      del = sum(sum(abs(av{n}(:,npt,:)-avp{n}(:,npt,:)).^2)); % errors
      mx = sum(sum(abs(av{n}(:,npt,:)).^2));
      if p.relerr && mx > p.tol
        del = del/mx;                            % normalise
      end
      e = del + e;                               % increment RMS error
  end                                            % end loop over averages
end                                              % end loop until time tmax
e = sqrt(e/(numel(p.averages)*p.points{1}(1)));  % get RMS mean error                                       
end                                              % end function