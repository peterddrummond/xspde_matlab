function [e,datn,in2,raw] = xcheck(checks,p)
%   [e,datn,in2,raw] = XCHECK(checks,r) checks different step-sizes.
%   It automatically runs xsim a total of 'checks' times, increasing
%   the initial r.steps by 2 after each run, to reduce the step-size by 2.
%   In each case it prints the maximum difference with any input 'compare'.
%   It also prints the statistical error-bar found at the maximum error.
%   The data returned is the data for xsim at the shortest stepsize.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

cl = fix(clock);                                    %%Start-up time
fprintf ('\nxCHECK, %d/%d/%d, time %d:%d:%d\n',cl); %%Print start-up time
p.checkmax = 1;
data = cell(checks);
t = cell(checks);
name = p.name;
[e,datn,in1,raw] = xsim(p);
data{1}=datn;
t{1}=e(6);
inputsteps = in1.steps(1);
for c = 2:checks
    p.steps  = inputsteps*2^(c-1);
    [e,datn,in2,raw] = xsim(p);
    data{c}=datn;
    t{c}=e(6);
end
xgraph(datn,in2);
dif = zeros(max(in2.averages),checks);
sg  = zeros(max(in2.averages),checks);
e1  = zeros(max(in2.averages),checks);
fprintf('\n======================================================+====\n');
fprintf(' Checks for %s',name);
for obs = in1.averages
  szd = size(data{1}{1}{obs});
  sz1 = prod(szd(1:end-1));
  sz2 = szd(end);
  if sz2 > 3
    if p.checkmax == 1
      fprintf('\n\n Max comparison errors in %s \n',in2.olabels{obs});
      fprintf(' \n Stepsize    Difference  Step error  Std dev.     Time');
      for c = 1:checks
        d1 = reshape(data{c}{1}{obs},[sz1,sz2]);
        [dif(obs,c),I] = max(abs((d1(:,1)-d1(:,4))));
        sg(obs,c)  = d1(I,3);
        e1(obs,c)  = d1(I,2);
        dt = in1.dt/(inputsteps*2^(c-1+in1.checks));
        fprintf('\n %.3e   %.3e   %.3e   %.3e    %.3e',...
        dt,dif(obs,c),e1(obs,c),sg(obs,c),t{c});
      end
    else
      fprintf('\n\n RMS comparison errors in %s \n',in2.olabels{obs});
      fprintf(' \n Stepsize    Difference  Step error  Std dev.     Time');
      for c = 1:checks
        d1 = reshape(data{c}{1}{obs},[sz1,sz2]);
        dif(obs,c) = sqrt(mean((d1(:,1)-d1(:,4)).^2));
        sg(obs,c)  = sqrt(mean(d1(:,3).^2));
        e1(obs,c)  = sqrt(mean(d1(:,2).^2));
        dt = in1.dt/(inputsteps*2^(c-1+in1.checks));
        fprintf('\n %.3e   %.3e   %.3e   %.3e    %.3e',...
        dt,dif(obs,c),e1(obs,c),sg(obs,c),t{c});
      end
    end
  end
end
fprintf('\n======================================================+====\n');
fprintf('\nxCHECK completed\n');
end