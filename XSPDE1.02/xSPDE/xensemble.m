function [data,raw] = xensemble (npar,latt) 
%   data = XENSEMBLE (npar,latt)  integrates and ensemble averages fields.
%   Input:  'npar' is total ensemble count,'latt' is  lattice description.
%   Output: 'data' is average data generated from stochastic integration.
%   and 'raw' is raw data generated from stochastic integration.
%   All xSPDE functions are licensed by Peter D. Drummond, (2015) - see License.txt 

sequence = length(latt);                          %%check length of input
r = latt{1};                                      %%first input
for seq=1:sequence                                %%loop over sequence
    data{seq} =  zeros(latt{seq}.d.data);         %%initialise data cells
end                                               %%end loop over sequence
raw = cell(r.ensembles(2),sequence,r.errorchecks);%%cell array for raw data
for ns = 1:r.ensembles(2)                         %%loop over ensembles
  nsp = (npar-1)*r.ensembles(2) + ns;             %%unique ensemble index 
  for  nc  = 1:r.errorchecks                      %%loop over errorchecks 
    rng(r.seed+nsp);                              %%Set unique random seed
    for seq = 1:sequence                          %%Loop over sequence
      if r.print                                  %%If print switch
        fprintf('Check %d, Ensemble %d\n',nc,nsp);%%indices
        if seq > 1                                %%If sequence
          fprintf('Sequence %d\n',seq);           %%print sequence number
        end;                                      %%end if sequence
      end;                                        %%end if print switch
      r = latt{seq};                              %%calculate grid 
      r.propagator = r.propfactor (nc,r);         %%calculate propagator   
      w = r.randomgen(r);                         %%generate random fields
      if seq == 1                                 %%check if sequence = 1
          a = r.initial(w,r);                     %%initialize fields
      else                                        %%if sequence > 1
          a = r.transfer(w,r,a,r0);               %%interface fields
      end                                         %%end check sequence =1
      r0 = r;                                     %%Old grid stored
      a = reshape(a,r.d.a);                       %%reshape fields
      r.t=r.origin(1);                            %%initial time
      [a,o,raw{ns,seq,nc}] = xpath(a,nc,r);       %%simulate path
      if nc == 2 || r.errorchecks == 1            %%if fine check
        data{seq}(1,:,:,:) = data{seq}(1,:,:,:) + o/r.n.ensemble;   %%fine means           
        data{seq}(3,:,:,:) = data{seq}(3,:,:,:) + o.^2/r.n.ensemble;%%fine var                                    
      else                                        %%else coarse calc
        data{seq}(2,:,:,:) = data{seq}(2,:,:,:) + o/r.n.ensemble;   %%coarse mean
      end                                         %%end if errorchecks
    end;                                          %%end sequence loop
  end                                             %%end errorchecks loop 
end;                                              %%end ensemble loop
end                                               %%end function
