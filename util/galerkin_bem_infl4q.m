function [He,Ge]=galerkin_bem_infl4q(ex1,ey1,ez1,ex,ey,ez,ep,n)
% [He,Ge]=bem_infl4q(coord,ex,ey,ez,ep,n)
% [He,Ge]=bem_infl4q(coord,ex,ey,ez,ep)
%----------------------------------------------------------------------
% PURPOSE
% Compute the element influence matrices He and Ge for a three
% dimensional four-node quadrilateral acoustic element.
%
% INPUT: coord=[x y z] coordinates of the influenced node
%
% ex=[x1 x2 x3 x4]
% ey=[y1 y2 y3 y4]
% ez=[z1 z2 z3 z4] node coordinates for the influencing
% element.
%
% ep = [w c rho] problem properties
% w: angular frequency
% c: speed of sound in acoustic medium
% rho: density of acoustic medium
%
% n=[value] normal direction
% value=1 default
% -1 reverse
%
% OUTPUT: He, Ge: Element influence matrices
%----------------------------------------------------------------------
rev=1;
if nargin==6
rev=n;
end
k=ep(1)/ep(2);
%****Gauss points****
ga=0.577350269189626;
xi=[-ga; ga; ga; -ga]; eta=[-ga; -ga; ga; ga];
N(:,1)=1/4*(xi-1).*(eta-1); N(:,2)=-1/4*(xi+1).*(eta-1);
N(:,3)=1/4*(xi+1).*(eta+1); N(:,4)=-1/4*(xi-1).*(eta+1);
xg=N*ex'; 
yg=N*ey'; 
zg=N*ez';
xf=N*ex1'; 
yf=N*ey1'; 
zf=N*ez1';

%****Element Area****

dNr(1:2:7,1)=-(1-eta)/4; dNr(1:2:7,2)= (1-eta)/4;
dNr(1:2:7,3)= (1+eta)/4; dNr(1:2:7,4)=-(1+eta)/4;
dNr(2:2:8,1)=-(1-xi)/4; dNr(2:2:8,2)=-(1+xi)/4;
dNr(2:2:8,3)= (1+xi)/4; dNr(2:2:8,4)= (1-xi)/4;
JTxy=dNr*[ex;ey]'; JTyz=dNr*[ey;ez]'; JTzx=dNr*[ez;ex]';
detJxy=[det(JTxy(1:2,:));det(JTxy(3:4,:));det(JTxy(5:6,:));
det(JTxy(7:8,:))];
detJyz=[det(JTyz(1:2,:));det(JTyz(3:4,:));det(JTyz(5:6,:));
det(JTyz(7:8,:))];
detJzx=[det(JTzx(1:2,:));det(JTzx(3:4,:));det(JTzx(5:6,:));
det(JTzx(7:8,:))];
A=[sqrt(detJxy.^2+detJyz.^2+detJzx.^2)];

JTxy=dNr*[ex1;ey1]'; JTyz=dNr*[ey1;ez1]'; JTzx=dNr*[ez1;ex1]';
detJxy=[det(JTxy(1:2,:));det(JTxy(3:4,:));det(JTxy(5:6,:));
det(JTxy(7:8,:))];
detJyz=[det(JTyz(1:2,:));det(JTyz(3:4,:));det(JTyz(5:6,:));
det(JTyz(7:8,:))];
detJzx=[det(JTzx(1:2,:));det(JTzx(3:4,:));det(JTzx(5:6,:));
det(JTzx(7:8,:))];
A1=[sqrt(detJxy.^2+detJyz.^2+detJzx.^2)];


%****Influence Vectors****
g = zeros(4,4);
h = zeros(4,4);
for i = 1:4
    for j = 1:4
      xdis=xg(j)-xf(i); ydis=yg(j)-yf(i); zdis=zg(j)-zf(i);
      dis=sqrt(xdis.^2+ydis.^2+zdis.^2);
      g(i,j)=1i*ep(3)*ep(1)*exp(-1i*k*dis)./(4*pi*dis)*A(j)*A1(i);
      h1=-xdis.*exp(-1i*k*dis)./(4*pi*dis.^2).*(1i*k+1./dis);
      h2=-ydis.*exp(-1i*k*dis)./(4*pi*dis.^2).*(1i*k+1./dis);
      h3=-zdis.*exp(-1i*k*dis)./(4*pi*dis.^2).*(1i*k+1./dis);
      a=[ex(2)-ex(1) ey(2)-ey(1) ez(2)-ez(1)]; b=[ex(4)-ex(1) ey(4)-ey(1) ez(4)-ez(1)]; n=[a(2)*b(3)-a(3)*b(2);a(3)*b(1)-a(1)*b(3); a(1)*b(2)-a(2)*b(1)];
      n=rev*n/sqrt(n'*n);
      h(i,j)=[h1 h2 h3]*n*A(j)*A1(i);
    end
end
Ge = N'*g*N;
He=-N'*h*N;
%-----------------------------------end--------------------------------