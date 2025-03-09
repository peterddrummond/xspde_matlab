function [a,av,raw] = xpath(a,nc,p)
%   [a,av,raw] = XPATH (a,nc,p)  integrates stochastic paths
%   Input:  a is a cell array of initial parallel trajectory field arrays
%           nc is the check index: nc = 1 for fine, nc = 2 for coarse
%           If nc = 1, and p.checks = 1, the checks loop is taken twice
%           p is the structure of parameters used
%   Output: a is a cell array of final trajectory field arrays
%           av is the average data generated during integration
%           raw is the raw data of all fields, if required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Called by xensemble
%   Calls method, define, xdata, breed, creshape, czeros
%   Licensed by Peter D. Drummond, (2024) - see License.txt, XSPDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

av  = czeros(p.d.av);                            % initialize average
raw = czeros(p.d.raw);                           % initialize raw cells
aft = czeros(p.d.ft);                            % initialize transforms
ast = czeros(p.d.ca);                            % initialize a-storage
p.breedw = 0;                                    % Initial breed count
inf = p.ind{1:p.nfields};
w = [];                                          % Set noise to void
n1 = 0;                                          % Initialise point count
f = p.fieldcells;                                % number of field cells
p.specpts =  p.plotpts+floor((p.steps(1)*p.points{1}(1))/2)+1;

%     LOOP IN TIME-POINTS, STEPS AND CHECKS, AND GENERATE NOISE TERMS

for n = 1:p.points{1}(1)                         % loop over time points 
 for step = 1:p.steps(1)                         % loop over steps
  n1 = 1+n1;                                     % Increment total count
  if n > 1 || step >1 || p.points{1}(1) == 1     % If not first point
    for checks = 0:p.checks(1)                   % loop if there are checks                      
      w = p.noisegen(w,nc,p.nscale,p);           % Compute noise terms
      if (nc == 1 || checks == 1)                % If not coarse t-check
        for ac = 1:p.auxcells                    % loop over auxiliaries
          a{ac+f} = p.define{ac}(a{1:f},w{:},p); % Get temporary auxiliary
        end                                      % end loop on auxiliary
        if p.jump{1} == 1
          a(1:f) = p.jmethod(a(1:f),w,p);        % jump integrate
        end
        a(1:f) = p.method(a(1:f),w,p);           % integrate field cells
        for ac = 1:p.auxcells                    % loop over auxiliaries
          a{ac+f} = 0.5*(a{ac+f}+p.define{ac}(a{1:f},w{:},p)); % Auxiliary
        end                                      % end loop on auxiliary
        if p.spectrum                            % If spectrum needed
          for l = 1:p.totcells                   % loop on total cells
            if checks == 1 && nc == 1
              aft{l}(inf,n1,p.ind{:}) = (reshape(a{l},p.d.a1{l})...
              + aft{l}(inf,n1,p.ind{:}))/2;
            else
              aft{l}(inf,n1,p.ind{:}) = reshape(a{l},p.d.a1{l});
            end
          end                                    % end loop on total cells
        end                                      % End if spectrum
        w = [];                                  % Clear used noise cells        
        p.t = p.t + p.dtr;                       % Update current time
      end                                        % end if high res or last     
    end                                          % end check loop
    if p.thresholdw > 0                          % Check if renormalizing
      [a{1:f},p.breedw] = p.breed(a{1:f},p);     % Renormalize and breed
    end                                          % End check for threshold
  else                                           % else if first point
  if p.spectrum                                  % If spectrum needed
    for l = 1:p.totcells                         % loop on total cells
      aft{l}(p.ind{1:p.nfields},n1,p.ind{:}) = reshape(a{l},p.d.a1{l});
    end                                          % end loop on total cells
  end                                            % end if spectrum
  end                                            % end if not first point
  
%     STORE RAW FIELDS, SPECTRAL FIELDS, TIME-DOMAIN OBSERVE DATA

  if step == 1                                   % check if the first step 
  for l = 1:p.totcells                           % loop on total cells  
    if p.rawdata                                 % check if raw data
      raw{l}(p.ind{1:p.nfields},n,p.ind{:}) = reshape(a{l},p.d.a1{l});% raw 
    end                                          % end if raw data needed
  end                                            % end cell loop
  for o = p.averages                             % loop over the averages
    if ~p.transforms{o}(1)                       % If not transformed
      av{o}(:,n,:) = xdata(a,o,p);               % average observables
    end                                          % End check if not
  end                                            % end loop on averages
  end                                            % end check if first 
 end                                             % end step loop
end                                              % end time loop
 
%     FOURIER TRANSFORM IN TIME AND STORE SPECTRAL OBSERVE DATA

if p.spectrum                                    % if spectrum is needed
  for l = 1:p.totcells                           % loop on total cells
    aft{l} = ifft(aft{l},p.points{1}(1)*p.steps,2)*p.fsc;%% Transform 
    aft{l} = fftshift(aft{l},2);                 % Use plot coordinates
    aft{l} = aft{l}(:,p.specpts,:); 
  end                                            % End loop on total cells  
  for n = 1:p.points{1}(1)                       % loop over points
    for l = 1:p.totcells                         % loop over total cells  
      ast{l} = reshape(aft{l}(:,n,:),p.d.ca{l}); % Reshape Fourier data
    end                                          % End loop on total cells  
    for o=p.averages                             % loop over the averages
      if p.transforms{o}(1)                      % Check if transformed
        av{o}(:,n,:) = xdata(ast,o,p);           % get averages if needed 
      end                                        % End check transformed
    end                                          % end loop over averages
  end                                            % end loop over points
end                                              % end if spectrum
end                                              % End function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XPATH FUNCTION