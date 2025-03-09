function [a,varargout]  =  RKWP21(a,w,p)        
%   a = RKWP21(a,w,p)  propagates an Ito step with 2nd-order weak Runge-Kutta. 
%   Treats a single field and noise, see Kloeden and Platen.
%   Uses derivA and derivB functions for the A and B functions.
%   The derivA and derivB functions are defined in the Ito calculus
%   where the noise term w is delta-correlated, so that <w^2> = 1/dt:
%
%   da/dt = derivA(a,p) + derivB(a,p)*w
%
%   Input: field a, noise w, parameters p. 
%   Output: new field a. 
%   Licensed by Run Yan Teh and Peter D. Drummond, (2024) - see License 

A   = p.derivA(a{1},p);
B   = p.derivB(a{1},p);
dt  = p.dtr; sdt = sqrt(dt);
ay  = a{1} + A*dt;
ap  = ay + B*sdt;
am  = ay - B*sdt;
ay  = ay + B.*w{1}*dt;
Bwp = p.derivB(ap,p);
Bwm = p.derivB(am,p);
a{1}   = a{1} + 0.5*(A + p.derivA(ay,p))*dt;
a{1}   = a{1} + 0.25*(Bwp + Bwm + 2*B).*w{1}*dt;
a{1}   = a{1} + 0.25*(Bwp - Bwm).*(w{1}.^2*dt - 1)*sdt;
varargout = {2,2,0,0,0}; 
%% St. order = 2, det. order = 2, ipsteps = 0, vect = 0(Y), cell  = 0(Y)
end