function z = xnoisefb(p)
%   z = XNOISEFB(p) stores total random noise 
%   Used for solving FBSDEs by iteration
%   Input:  p is the structure of parameters used
%   Output: z is a nested cell array of noise trajectories
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Called by xpathfb
%   Needs: steps, checks, points, noisecells, p.d.noise, nscale, noisefactor
%   Licensed by Peter D. Drummond, (2025) - see License.txt, XSPDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nt=p.steps(1)*(p.checks(1)+1)*(p.points{1}(1)-1);% Number of time points
z = cell(1,nt);                                  % Gaussian noise cells
for np = 1:nt                                    % loop over time points
  for c = 1:p.noisecells(1)                      % Loop over noise cells
    w = randn(p.d.noise{c});                     % Gaussian noises
    w = w*sqrt(p.nscale*p.noisefactor{c});       % Gaussian scale factor
    z{np}{c} = w;                                % store the  noise
  end                                            % end loop over cells
end                                              % end loop over time
end                                              % end function