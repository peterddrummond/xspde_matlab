function [data,raw] = xensemble (npar,l) 
%   data = XENSEMBLE (npar,latt)  integrates and ensemble averages fields.
%   Input:  'npar' is the ensemble number,'l' is the simulation parameters.
%   Output: 'data' is average data generated from stochastic integration.
%           'raw' is raw data generated from stochastic integration.
%   First dimension is the field index, last dimension is the ensemble
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License


sequence = length(l);                             %%check length of input
data = cell(sequence);                            %%allocate data sequence
for seq=1:sequence                                %%loop over sequence
    for n=1:l{seq}.functions
        data{seq}{n} = zeros(l{seq}.d.data{n});   %%initialise data cells
    end
end                                               %%end loop over sequence
serial = l{1}.ensembles(2);                       %%serial ensemble number
raw = cell(sequence,l{1}.errorchecks,serial);     %%cell array for raw data

      %%Loop over all the serial ensembles

for ns = 1:serial                                 %%loop over ensembles
  nsp = (npar-1)*serial + ns;                     %%unique ensemble index  
  for  nc  = 1:l{1}.errorchecks                   %%loop over errorchecks
    if l{1}.octave
            randn('state',l{1}.seed+nsp)
    else
            rng(l{1}.seed+nsp);
    end                                           %%Set unique random seed
    if l{1}.print                                 %%If print switch
        fprintf('Check %d, Ensemble %d, Time %f\n',nc,nsp,toc());%%indices
    end                                           %%end if print switch
    
      %%Loop over the stochastic integration sequence   
    
    for seq = 1:sequence                          %%Loop over sequence
      r=l{seq};
      if r.print > 1                              %%If verbose print
          fprintf('Sequence %d\n',seq);           %%print sequence indices
      end                                         %%end if verbose  print
      r.dtr=r.dt/nc;                              %%reduced step-size
      r.propagator = r.propfactor(nc,r);          %%get propagator  
      w = r.randomgen(r);                         %%generate random fields
      if seq == 1                                 %%check if sequence = 1
          if r.setboundaries                      %%initialize boundaries
              r.t  = r.origin(1)-1;               %%Set time < origin
              for dir=2:r.dimension
                  r.boundval{dir} = r.boundfun(w,dir,r); %%Store boundaries
              end
          end
          a = r.initial(w,r);                     %%initialize fields
          if length(a) == 1
              a = a*ones(r.d.a);
          end
      else                                        %%if sequence > 1
          [a,r] = r.transfer(w,r,a(1:r.fields,:),l{seq-1}); %%interface fields
      end                                         %%end check sequence =1
      a = reshape(a,r.d.a);                       %%reshape fields
      
      r.t=r.origin(1);                            %%initial time
      
      %%Store the averages for the stochastic field path
      
      if r.transformw
           [a,av,raw{seq,nc,ns}] = xpathw(a,nc,r);%%simulate path
      else
           [a,av,raw{seq,nc,ns}] = xpath(a,nc,r); %%simulate path
      end 
      for n = 1:r.averages
          av{n} = reshape(av{n},r.d.av{n});       %%reshape average data
      end                                         %%end averages loop
      for n = 1:r.functions
        f = r.function{n}(av,r);                  %%graph data
        f = reshape(f,r.d.sf{n});                 %%reshape graph data  
        if nc == 2 || r.errorchecks == 1          %%if fine check
         data{seq}{n}(:,1) = data{seq}{n}(:,1) + f/r.ncopies;          
         data{seq}{n}(:,3) = data{seq}{n}(:,3) + f.*conj(f)/r.ncopies;                                 
        else                                      %%else coarse calc
         data{seq}{n}(:,2) = data{seq}{n}(:,2) + f/r.ncopies;
        end                                       %%end if errorchecks
      end                                         %%end functions loop
    end                                           %%end sequence loop
  end                                             %%end errorchecks loop 
end                                               %%end ensemble loop
end                                               %%end function