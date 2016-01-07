function [data,raw] = xensemble (npar,l) 
%   data = XENSEMBLE (npar,latt)  integrates and ensemble averages fields.
%   Input:  'npar' is total ensemble count,'l' is  lattice description.
%   Output: 'data' is average data generated from stochastic integration.
%   and 'raw' is raw data generated from stochastic integration.
%   xSPDE functions are licensed by Peter D. Drummond, (2015) - see License

sequence = length(l);                             %%check length of input
data = cell(sequence);                            %%allocate data sequence
for seq=1:sequence                                %%loop over sequence
    data{seq} =  zeros(l{seq}.d.data);            %%initialise data cells
end                                               %%end loop over sequence
serial = l{1}.ensembles(2);                       %%serial ensemble number
raw = cell(serial,sequence,l{1}.errorchecks);     %%cell array for raw data
for ns = 1:serial                                 %%loop over ensembles
  nsp = (npar-1)*serial + ns;                     %%unique ensemble index 
  for  nc  = 1:l{1}.errorchecks                   %%loop over errorchecks 
    rng(l{1}.seed+nsp);                           %%Set unique random seed
    if l{1}.print                                 %%If print switch
        fprintf('Check %d, Ensemble %d\n',nc,nsp);%%print indices
    end;                                          %%end if print switch
    for seq = 1:sequence                          %%Loop over sequence
      if l{seq}.print > 1                         %%If verbose print
          fprintf('Sequence %d\n',seq);           %%print sequence indices
      end;                                        %%end if verbose  print
      l{seq}.propagator = l{seq}.propfactor(nc,l{seq});%%get propagator   
      w = l{seq}.randomgen(l{seq});               %%generate random fields
      if seq == 1                                 %%check if sequence = 1
          a = l{seq}.initial(w,l{seq});           %%initialize fields
      else                                        %%if sequence > 1
          a = l{seq}.transfer(w,l{seq},a,l{seq-1});%%interface fields
      end                                         %%end check sequence =1
      a = reshape(a,l{seq}.d.a);                  %%reshape fields
      l{seq}.t=l{seq}.origin(1);                  %%initial time
      [a,o,raw{ns,seq,nc}] = xpath(a,nc,l{seq});  %%simulate path
      if nc == 2 || l{seq}.errorchecks == 1       %%if fine check
        data{seq}(1,:,:,:) = data{seq}(1,:,:,:) + o/l{seq}.ncopies; %%means           
        data{seq}(3,:,:,:) = data{seq}(3,:,:,:) + o.^2/l{seq}.ncopies;%%var                                    
      else                                        %%else coarse calc
        data{seq}(2,:,:,:) = data{seq}(2,:,:,:) + o/l{seq}.ncopies;%%coarse
      end                                         %%end if errorchecks
    end;                                          %%end sequence loop
  end                                             %%end errorchecks loop 
end;                                              %%end ensemble loop
end                                               %%end function

%   Version 1.03 - improved ensemble and sequence progress printout
