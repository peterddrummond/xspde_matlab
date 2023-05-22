function [e]  =  Boundaries2()
%   [e]  =  BOUNDARIES is a test of xSPDE boundary types.
%   Solves a (one+one) dimensional partial differential equation
%
%   Dataset (I) uses a trigonometric transform linear propagater term
%   Dataset (II) uses finite difference methods
%   
%   There are five tests in each dataset, for a five-component field
%
%   (1) Using Dirichlet or zero field boundary conditions
%   (2) Using Robin or zero gradient boundary conditions
%   (3) Using Dirichlet-Robin boundary conditions
%   (4) Using Robin-Dirichlet boundary conditions
%   (5) Using periodic boundary conditions
%
%   All results are compared to exact solutions for a heat equation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.dimensions     =  2;
p.points         =  [51,51];
p.order          =  0;
p.verbose        =  1;
p.fields         =  5;
p.ranges         =  [4,pi];
p.origins        =  [0,0];
p.initial        =  @heat_in; 
p.observe        =  {@(a,p) a(1,:),@(a,p) a(2,:),@(a,p) a(3,:), ...
                     @(a,p) a(4,:),@(a,p) a(5,:)};
p.compare        =  {@heat_1,@heat_2,@heat_3,@heat_4,@heat_5};
p.diffplot       =  {1,1,1,1,1};
p.olabels        =  {'a, DD','a, NN','a, DN','a, ND','a, PP'};
p.name           =  'Boundaries2: 1+1D spectral';
p.boundaries{2}  =  [1,1;-1,-1;1,-1;-1,1;0,0];
p.linear         =  @(p) p.Dx.^2;
p1               =  p;
p1.linear        =  @(p) [];
p1.deriv         =  @(a,w,p) D2(a,2,p);
p1.steps         =  40;
p1.transfer      =  @(~,p,~,~) heat_in(0,p);
p1.name          =  'Boundaries2: 1+1D finite diffs';
e                =  xspde({p,p1});                      
end                                                       

function a  =  heat_in(~,p)
  a(1,:) = 4*sin(p.x)+sin(2*p.x); 
  a(2,:) = 5+4*cos(p.x)+cos(2*p.x);
  a(3,:) = 4*sin(p.x/2)+sin(3*p.x/2); 
  a(4,:) = 4*cos(p.x/2)+cos(3*p.x/2);
  a(5,:) = 2+cos(2*p.x/1.02)+sin(4*p.x/1.02);
end

function o  =  heat_1(p)
  o = 4*sin(p.x).*exp(-p.t)+sin(2*p.x).*exp(-4*p.t); 
end
function o  =  heat_2(p)
  o = 5+4*cos(p.x).*exp(-p.t)+cos(2*p.x).*exp(-4*p.t);  
end
function o  =  heat_3(p)
  o = 4*sin(p.x/2).*exp(-p.t/4)+sin(3*p.x/2).*exp(-9*p.t/4);  
end
function o  =  heat_4(p) 
  o = 4*cos(p.x/2).*exp(-p.t/4)+cos(3*p.x/2).*exp(-9*p.t/4);  
end
function o  =  heat_5(p) 
  o = 2+cos(2*p.x/1.02).*exp(-4*p.t/1.02^2)+...
      sin(4*p.x/1.02).*exp(-16*p.t/1.02^2);  
end
