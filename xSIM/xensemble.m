function [data,raw] = xensemble (npar,input) 
%   [data,raw] = XENSEMBLE (npar,input)  integrates stochastic ensembles
%   'npar' is the ensemble number,'input' is the parameter cell array.
%   Output: 'data' is average data generated from stochastic integration.
%           'raw' is raw data of stochastic plus auxiliary fields
%   Called by xsim
%   Licensed by Peter D. Drummond, (2023) - see License.txt, XSDE manual 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

sequence = length(input);                        % check length of input
data   = cell(sequence);                         % allocate data sequence
p1   = input{1};                                 % first parameter struct
raw = cell(sequence,p1.nchk,p1.ensembles(2));    % initialize raw cells
D   = cell(sequence);
ne = 1./(p1.ensembles(2)*p1.ensembles(3)) - 1;   % ensemble scaling factor
nch1 = max(3,1+p1.nchk);
for s = 1:sequence
  for o = 1:input{s}.outputs
    data{s}{o}  = zeros([input{s}.d.data{o},nch1]);
  end
end
for j=1:p1.ensembles(2)                          % loop over ensembles
  nsp = p1.seed+npar+(j-1)*p1.ensembles(3);      % Unique random seed 
  for  nc  = 1:p1.nchk                           % loop over errorchecks
    if p1.octave                                 % check if octave
      randn('state',p1.seed+nsp)                 %#ok<RAND> %Set the seed
    else                                         %else matlab
      rng(nsp,p1.rng);                           % Set unique random seed  
    end                                          % end check if octave
    
%     LOOP OVER THE INPUT SEQUENCE                                                  
 
    for s = 1:sequence                           % loop over sequence
      p   =    input{s};                         % get the parameters
      if p.verbose > 0                           % if verbose printing
        fprintf('Ensemble %d.%d, check %d, seq %d\n',npar,j,nc,s);
      end                                        % end if verbose printing
      if p.checks(1) && nc ~= 2
        p.dtr =  p.dt/(2*p.steps(1));            % reduced step-size in t
      else
        p.dtr =  p.dt/(p.steps(1));              % step-size in t
      end
      a     = cell(1,p.fieldcells);              % initialize field cells
      p.propagator = cell(1,p.fieldcells);       % initialize propagator
      p.da = cell(1,p.fieldcells);               % initialize da cells
      for c = 1:p.fieldcells                     % loop over field cells
        if isequal (p.linear{c},0) || isequal (p.linear{c},[])
          p.propagator{c} = 1;
        else
          [p.propagator{c},p.da{c}] = p.propfactor(c,p); % get propagator
        end
        a0  =  zeros(p.d.ca{c});                 % initialize a0 to zeros
      end
      v = p.randomgen(p);
      for c = 1:p.fieldcells                     % loop over field cells    
        if s == 1                                % If first in sequence
          a{c} = p.initial{c}(v{:},p);           % initial  vectors
        else                                     % Else not first
          a{c} = p.transfer{c}(a1{:},v{:},p);    % transfer to next
        end                                      % End if first
        a{c} = xshape(a{c},c,0,p)+a0;            % reshape fields
        if p.dimensions > 1 && p.setboundaries   % If any boundaries
          p.t  = p.origins(1)-1;                 % Set time < origin
          for e = 2:p.dimensions                 % Loop over dimension
           p.boundval{c,e} = p.boundfun(a,c,e,p);% Store boundary
          end                                    % End loop over dimension
        end                                      % end initial boundaries
      end                                        % end loop over cells
      for ac = 1:p.auxcells                      % if auxiliary cells
        w = xnoise([],nc,p.nscale,p);            % add noise for defines
        a = [a,p.define{ac}(a{1:p.fieldcells},w{:},p)];%#ok<AGROW>
      end                                        % end if auxiliary cells
      p.t   =  p.origins(1);                     % set time origin
      [a1,av,raw{s,nc,j}] = p.path(a,nc,p);      % compute the path   
      
%     LOOP OVER THE OUTPUTS     
    
      for o=1:p.outputs                          % start outputs loop
        if p.d.data{o}(1) > 0                    % if data defined
          nes = 1+isequal(p.scatters{o},0)*ne;   % get scaling factor
          if isequal(p.scatters{o},0) || j < 2   % non-scatters or first
            D{s}{nc,o} = reshape(real(p.output{o}(av,p)),[p.d.data{o},1]); 
            if nc == p.nchk                      % wait till last check                      
              if p1.checks(1) &&  p.order > 0    % extrapolating in time? 
                Del =  (D{s}{1,o} - D{s}{2,o})/(2^(p.order)-1.);
                D{s}{2,o} = D{s}{1,o};
                D{s}{1,o} = D{s}{1,o}+Del;
              end                                % end if extrapolate
              for n1 = 1:p.nchk                  % loop over checks
                data{s}{o}(:,n1) = data{s}{o}(:,n1) + D{s}{n1,o}*nes;
              end
              data{s}{o}(:,nch1) = data{s}{o}(:,nch1)+D{s}{1,o}.^2*nes;
            end                                  % end if last check 
          end                                    % end if data needed
        end                                      % end if data defined
      end                                        % end functions loop
    end                                          % end sequence loop
  end                                            % end checks loop
end                                              % end ensemble loop
%%%%%%%%%%%%%%%%%%%%%%%  END XENSEMBLE FUNCTION