function [e]  =  Boundaries3()
%   [e]  =  BOUNDARIES3 is a test of xSPDE boundary types.
%   Solves a (one+two) dimensional stochastic partial differential equation
%   with a trigonometric transform linear propagater term
%   
%   There are six tests in each dataset, for a six-component field.
%
%   (1) Using Dirichlet in x, Robin in y
%   (2) Using Robin in x, Dirichlet in y
%   (3) Using Dirichlet-Robin boundary conditions in  x,y 
%   (4) Using Robin-Dirichlet boundary conditions in x,y
%   (5) Using periodic boundary conditions in dimensions x,y
%   (6) Using Dirichlet or zero field boundary conditions in x
%       and Dirichlet-Robin boundary conditions in y
%
%   All results are compared to exact solutions for a heat equation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.dimensions     =  3;
p.points         =  [21,21,21];
p.order          =  0;
p.fields         =  6;
p.noises         =  6;
p.ensembles      =  [8,1,8];
p.ranges         =  [4,pi,pi];
p.origins        =  [0,0,0];
p.initial        =  @heat_in; 
p.observe        =  {@(a,p) a(1,:),@(a,p) a(2,:),@(a,p) a(3,:),...
                     @(a,p) a(4,:),@(a,p) a(5,:),@(a,p) a(6,:)};
p.compare        =  {@heat_1,@heat_2,@heat_3,@heat_4,@heat_5,@heat_6};
p.diffplot       =  {1,1,1,1,1,1};
p.olabels        =  {'u, DD,RR','u, RR,DD','u, DR,DR','u, RD,RD',...
                     'u, PP,PP','u, DD,DR'};
p.name           =  'Boundaries3: 2+1D spectral';
p.boundaries{2}  =  [1,1;-1,-1;1,-1;-1,1;0,0;1, 1];
p.boundaries{3}  =  [-1,-1;1,1;1,-1;-1,1;0,0;1,-1];
p.deriv          =  @(a,w,p) 0.1*w;
p.linear         =  @(p) p.Dx.^2+p.Dy.^2;
e                =  xspde(p);               
end                                                       

function a  =  heat_in(~,p)
  a(1,:,:) = 4*sin(p.x).*cos(p.y)+sin(2*p.x).*cos(2*p.y); 
  a(2,:,:) = 4*cos(p.x).*sin(p.y)+cos(2*p.x).*sin(2*p.y);
  a(3,:,:) = 4*sin(p.x/2).*sin(p.y/2)+sin(3*p.x/2).*sin(3*p.y/2); 
  a(4,:,:) = 4*cos(p.x/2).*cos(p.y/2)+cos(3*p.x/2).*cos(3*p.y/2);
  a(5,:,:) = 2+cos(2*p.x/1.05).*cos(2*p.y/1.05)+sin(4*p.x/1.05).*sin(4*p.y/1.05);
  a(6,:,:) = 4*sin(p.x).*sin(p.y/2)+sin(2*p.x).*sin(3*p.y/2);
end

function o  =  heat_1(p)
  o = 4*sin(p.x).*cos(p.y).*exp(-2*p.t)+sin(2*p.x).*cos(2*p.y).*exp(-8*p.t); 
end
function o  =  heat_2(p)
  o = 4*cos(p.x).*sin(p.y).*exp(-2*p.t)+cos(2*p.x).*sin(2*p.y).*exp(-8*p.t);  
end
function o  =  heat_3(p)
  o = 4*sin(p.x/2).*sin(p.y/2).*exp(-p.t/2)+sin(3*p.x/2).*sin(3*p.y/2).*exp(-9*p.t/2);  
end
function o  =  heat_4(p) 
  o = 4*cos(p.x/2).*cos(p.y/2).*exp(-p.t/2)+cos(3*p.x/2).*cos(3*p.y/2).*exp(-9*p.t/2);  
end
function o  =  heat_5(p) 
  o = 2+cos(2*p.x/1.05).*cos(2*p.y/1.05).*exp(-8*p.t/1.05^2)+...
      sin(4*p.x/1.05).*sin(4*p.y/1.05).*exp(-32*p.t/1.05^2);  
end
function o  =  heat_6(p)
  o = 4*sin(p.x).*sin(p.y/2).*exp(-5/4*p.t)+sin(2*p.x).*sin(3*p.y/2).*exp(-25/4*p.t);
end
