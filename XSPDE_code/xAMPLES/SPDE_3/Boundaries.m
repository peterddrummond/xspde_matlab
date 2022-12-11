function [e]  =  Boundaries()
%   [e]  =  BOUNDARIES is a test of xSPDE boundary types.
%   Solves a (one+one) dimensional partial differential equation
%
%   Dataset (I) uses a trigonometric transform linear propagater term
%   Dataset (II) uses finite difference methods
%   
%   There are six tests in each dataset, for a five-component field
%
%   (1) Using Dirichlet or zero field boundary conditions in dimensions 1&2
%   (2) Using Robin or zero gradient boundary conditions in dimensions 1&2
%   (3) Using Dirichlet-Robin boundary conditions in dimensions 1&2 
%   (4) Using Robin-Dirichlet boundary conditions in dimensions 1&2
%   (5) Using periodic boundary conditions in dimensions 1&2
%   (6) Using Dirichlet or zero field boundary conditions in dimension 1
%       and Dirichlet-Robin boundary conditions in dimension 2
%
%   All results are compared to exact solutions for a heat equation
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p.dimensions     =  3;
p.points         =  [51,51,51];
p.order          =  0;
p.verbose        =  1;
p.fields         =  6;
p.ranges         =  [4,pi,pi];
p.origins        =  [0,0,0];
p.initial        =  @heat_in; 
p.observe        =  {@(a,p) a(1,:),@(a,p) a(2,:),@(a,p) a(3,:),...
                     @(a,p) a(4,:),@(a,p) a(5,:),@(a,p) a(6,:)};
p.compare        =  {@heat_1,@heat_2,@heat_3,@heat_4,@heat_5,@heat_6};
p.diffplot       =  {1,1,1,1,1,1};
p.olabels        =  {'u, DD,DD','u, RR,RR','u, DR,DR','u, RD,RD','u, PP,PP',...
                     'u, DD,DR'};
p.name           =  'Heat test, spectral';
p.boundaries{2}  =  [1,1;-1,-1;1,-1;-1,1;0,0;1, 1];
p.boundaries{3}  =  [1,1;-1,-1;1,-1;-1,1;0,0;1,-1];
p.linear         =  @(p) p.Dx.^2+p.Dy.^2;
p1               =  p;
p1.linear        =  @(p) [];
p1.deriv         =  @(a,w,p) D2(a,2,p)+D2(a,3,p);
p1.steps         =  200;
p1.transfer      =  @heat_tr; 
p1.name          =  'Heat test, finite diffs';
e                =  xspde({p,p1}); 
%e                =  xspde(p);                        

end                                                       

function a  =  heat_in(~,p)
  a(1,:,:) = 4*sin(p.x).*sin(p.y)+sin(2*p.x).*sin(2*p.y); 
  a(2,:,:) = 4*cos(p.x).*cos(p.y)+cos(2*p.x).*cos(2*p.y);
  a(3,:,:) = 4*sin(p.x/2).*sin(p.y/2)+sin(3*p.x/2).*sin(3*p.y/2); 
  a(4,:,:) = 4*cos(p.x/2).*cos(p.y/2)+cos(3*p.x/2).*cos(3*p.y/2);
  a(5,:,:) = 2+cos(2*p.x/1.02).*cos(2*p.y/1.02)+sin(4*p.x/1.02).*sin(4*p.y/1.02);
  a(6,:,:) = 4*sin(p.x).*sin(p.y/2)+sin(2*p.x).*sin(3*p.y/2);
end

function a  =  heat_tr(~,p,~,~)
  a  =  heat_in(0,p);
end


function o  =  heat_1(p)
  o = 4*sin(p.x).*sin(p.y).*exp(-2*p.t)+sin(2*p.x).*sin(2*p.y).*exp(-8*p.t); 
end
function o  =  heat_2(p)
  o = 4*cos(p.x).*cos(p.y).*exp(-2*p.t)+cos(2*p.x).*cos(2*p.y).*exp(-8*p.t);  
end
function o  =  heat_3(p)
  o = 4*sin(p.x/2).*sin(p.y/2).*exp(-p.t/2)+sin(3*p.x/2).*sin(3*p.y/2).*exp(-9*p.t/2);  
end
function o  =  heat_4(p) 
  o = 4*cos(p.x/2).*cos(p.y/2).*exp(-p.t/2)+cos(3*p.x/2).*cos(3*p.y/2).*exp(-9*p.t/2);  
end
function o  =  heat_5(p) 
  o = 2+cos(2*p.x/1.02).*cos(2*p.y/1.02).*exp(-8*p.t/1.02^2)+...
      sin(4*p.x/1.02).*sin(4*p.y/1.02).*exp(-32*p.t/1.02^2);  
end
function o  =  heat_6(p)
  o = 4*sin(p.x).*sin(p.y/2).*exp(-5/4*p.t)+sin(2*p.x).*sin(3*p.y/2).*exp(-25/4*p.t);
end
