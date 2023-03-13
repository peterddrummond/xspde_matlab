function [data,raw] = xensemble (npar,l) 
%   data = XENSEMBLE (npar,latt)  integrates and ensemble averages fields.
%   Input:  'npar' is the ensemble number,'l' is the parameter structure.
%   Output: 'data' is average data generated from stochastic integration.
%           'raw' is raw data of stochastic plus auxiliary fields
%   First dimension is the field index, last dimension is the ensemble
%   Called by: xsim
%   Needs: xpathfb,xpath,xpathw
%   Licensed by Peter D. Drummond, (2021) - see License.txt
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% INITIALIZE OUTPUT DATA CELLS
%
sequence = length(l);                             %%check length of input
data = cell(sequence);                            %%allocate data sequence
for seq=1:sequence                                %%loop over sequence
    for n=1:l{seq}.functions                      %%loop over outputs
        data{seq}{n} = zeros(l{seq}.d.data{n});   %%initialise data cells
    end                                           %%end loop over outputs
end                                               %%end loop over sequence
serial = l{1}.ensembles(2);                       %%serial ensemble number
parallel = l{1}.ensembles(3);                     %%serial ensemble number
raw = cell(sequence,l{1}.errorchecks,serial);     %%cell array for raw data

      %%Loop over all the serial ensembles

for ns = 1:serial                                 %%loop over ensembles
  nsp = (npar-1)*serial + ns;                     %%unique ensemble index
  p = l{1};
  for  nc  = 1:p.errorchecks                      %%loop over errorchecks
    if p.octave                                   %%check if octave
            randn('state',l{1}.seed+nsp)          %#ok<RAND> %%Set the seed
    else                                          %%else matlab
            rng(l{1}.seed+nsp);                   %%Set unique random seed
    end                                           %%End check if octave
    xpr(1,p,'Check:%d, Serial:%d/%d, Parallel:%d/%d, Time:%f\n',...
        nc,ns,serial,npar,parallel,toc());        %%indices
    
      %%Loop over the stochastic integration sequence   
    
    for seq = 1:sequence                          %%Loop over sequence
      p=l{seq};                                   %%get parameter struct
      xpr(2,p,'Sequence %d\n',seq);               %%print sequence indices
      p.dtr=p.dt/nc;                              %%reduced step-size
      p.propagator = p.propfactor(nc,p);          %%get propagator
      w = p.randomgen(p);                         %%generate random fields
      p.t = p.origin(1);                          %%initial time
      if p.fields > p.fieldsf                     %%check if fb case
          [a,av,raw{seq,nc,ns}] = xpathfb(w,nc,p);%%simulate fb path
      else                                        %%else standard case
        if seq == 1                               %%check if sequence = 1
          if p.setboundaries                      %%initialize boundaries
              p.t  = p.origin(1)-1;               %%Set time < origin
              for dir=2:p.dimension               %%Loop over dimension
                  p.boundval{dir} = p.boundfun(w,dir,p); %%Store boundaries
              end                                 %%End loop over dimension
          end                                     %%end initial boundaries
          p.t  = p.origin(1);                     %%Set time = origin
          a = p.initial(w,p);                     %%initialize fields
          a = xshape(a,0,p)+zeros(p.d.a);         %%reshape fields
        else                                      %%if sequence > 1
          a = p.transfer(w,p,a,l{seq-1});
          a = xshape(a,0,p)+zeros(p.d.a);         %%reshape fields
        end                                       %%end check sequence =1
        if p.transformw                           %%If spectrum needed
           [a,av,raw{seq,nc,ns}] = xpathw(a,nc,p);%%simulate spectral path
        else                                      %%else no spectrum needed
           [a,av,raw{seq,nc,ns}] = xpath(a,nc,p); %%simulate normal path    
        end                                       %%End if spectrum needed
      end                                         %%End check if fb case
      
      %%Reshape the averages for the stochastic field path
           
      for n = 1:p.averages                        %%loop over averages
          av{n} = reshape(av{n},p.d.av{n});       %%reshape average data
      end                                         %%end loop over averages
      
      %%Store the functions and comparisons for the stochastic field path
      
      p.indext = 1;                               %%switch for time index
      for n = 1:p.functions                       %%loop over outputs
        f = p.function{n}(av,p);                  %%get graphics data 
        f = reshape(f,p.d.sf{n});                 %%reshape graph data
        if  p.scatters{n} > 1                     %%if scatters
          if nc == 2 || p.errorchecks == 1        %%if fine check
            data{seq}{n}(:,1) = f;          
            data{seq}{n}(:,3) = f.*conj(f);                                 
          else                                    %%else coarse calc
            data{seq}{n}(:,2) = f;
          end                                     %%end if errorchecks
        else
          if nc == 2 || p.errorchecks == 1        %%if fine check
            data{seq}{n}(:,1) = data{seq}{n}(:,1) + f/p.ncopies;          
            data{seq}{n}(:,3) = data{seq}{n}(:,3) + f.*conj(f)/p.ncopies;                                 
          else                                    %%else coarse calc
            data{seq}{n}(:,2) = data{seq}{n}(:,2) + f/p.ncopies;
          end                                     %%end if scatters
        end                                       %%end if errorchecks
      end                                         %%end functions loop
    end                                           %%end sequence loop
  end                                             %%end errorchecks loop 
end                                               %%end ensemble loop
end                                               %%end function