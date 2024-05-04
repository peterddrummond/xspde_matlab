function [a,av,raw] = xpath(a,nc,p)
%   [a,av,raw] = XPATH (a,nc,p)  integrates stochastic paths
%   Input:  a is a cell array of initial parallel trajectory field arrays
%           nc is the check index: nc = 1 for coarse, nc = 2 for fine
%           p is the structure of parameters used
%   Output: a is a cell array of final trajectory field arrays
%           av is the average data generated during integration
%           raw is the raw data of all fields, if required
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Called by xensemble
%   Calls method, define, xdata, breed, creshape, czeros
%   Licensed by Peter D. Drummond, (2024) - see License.txt, XSPDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

av  = czeros(p.d.av);                            % initialize averages
raw = czeros(p.d.raw);                           % initialize raw cells
aft = czeros(p.d.ft);                            % initialize transforms
ast = czeros(p.d.ca);                            % initialize a-storage
p.breedw = 0;                                    % Initial breed count
x = 1/2;                                         % Auxiliary constant
if p.checks(1) && nc ~= 2
    x = 1/4;
end
w = [];                                          % Set noise to void
f = p.fieldcells;                                % number of field cells
ax1 = czeros(p.d.ca(1+f:end));                   % initial auxiliary temp

%     LOOP IN TIME-POINTS, STEPS AND CHECKS, AND GENERATE NOISE TERMS

for n = 1:p.points{1}(1)                         % loop over time points
 if n > 1  ||  p.points{1}(1) == 1               % If not first point
  for step = 1:p.steps(1)                        % loop over steps
   ax = czeros(p.d.ca(1+f:end));                 % initialize auxiliary
   for checks = 0:p.checks(1)                    % loop on checks                      
    w = p.noisegen(w,nc,p.nscale,p);             % Compute noise terms
    %if nc > 1 || checks == p.checks(1)          % If hi res or last check
    if nc ~= 2 || ~p.checks(1) || checks ~= 0    % If not coarse t-check
      for l = 1:p.auxcells                       % loop over auxiliaries
        ax1{l} = p.define{l}(a{l},w{l},p);       % Get temporary auxiliary
        ax{l} = ax{l}+x*ax1{l};                  % Add temp to auxiliary
      end                                        % end loop on auxiliary
      if p.jump{1} == 1
          a(1:f) = p.jmethod(a(1:f),w,p);        % jump integrate
      end
      a(1:f) = p.method(a(1:f),w,p);             % integrate cells
      for l = 1:p.auxcells                       % loop on auxiliaries   
        ax1{l} = p.define{l}(a{l},w{l},p);       % Get temporary auxiliary
        ax{l} = ax{l}+x*ax1{l};                  % Add temp to auxiliary
      end                                        % end loop on auxiliaries   
      w = [];                                    % Clear used noise cells        
      p.t = p.t + p.dtr;                         % Update current time
    end                                          % end if high res or last     
   end                                           % end check loop  
  end                                            % end step loop

%     STORE RAW FIELDS, SPECTRAL FIELDS, TIME-DOMAIN OBSERVE DATA

  a(f+1:end) = ax;                               % Store auxiliary 
  if p.thresholdw > 0                            % Check if renormalizing
    [a{1:f},p.breedw] = p.breed(a{1:f},p);       % Renormalize and breed
  end                                            % End check for threshold
 end                                             % end if np > 1
 for l = 1:p.totcells                            % loop on total cells
    if p.spectrum                                % If spectrum needed
      aft{l}(:,n,:) = reshape(a{l},p.d.a1{l});   % reshape and store data 
    end                                          % End if spectrum
    if p.rawdata                                 % check if raw data 
      raw{l}(:,n,:) = reshape(a{l},p.d.a1{l});   % store raw data
   end                                           % end if raw data needed
 end                                             % end cell loop
 for o = 1:p.averages                            % loop over the averages
    if ~p.transforms{o}(1)                       % If not transformed
      av{o}(:,n,:) = xdata(a,o,p);               % get average observables
    end                                          % End check if not
 end                                             % end loop on averages
end                                              % end time loop
 
%     FOURIER TRANSFORM IN TIME AND STORE SPECTRAL OBSERVE DATA

if p.spectrum                                    % if spectrum is needed
  for l = 1:p.totcells                           % loop on total cells
    aft{l} = ifft(aft{l},p.points{1}(1),2)*p.fsc;%% Transform in time
    aft{l} = fftshift(aft{l},2);                 % Use plot coordinates
  end                                            % End loop on total cells  
  for n = 1:p.points{1}(1)                       % loop over points
    for l = 1:p.totcells                         % loop over total cells  
      ast{l} = reshape(aft{l}(:,n,:),p.d.ca{l}); % Reshape Fourier data
    end                                          % End loop on total cells  
    for o=1:p.averages                           % loop over the averages
      if p.transforms{o}(1)                      % Check if transformed
        av{o}(:,n,:) = xdata(ast,o,p);           % get averages if needed 
      end                                        % End check transformed
    end                                          % end loop over averages
  end                                            % end loop over points
end                                              % end if spectrum
end                                              % End function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  END XPATH FUNCTION