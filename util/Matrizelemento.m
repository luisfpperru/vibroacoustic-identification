function [N,A,Area]= Matrizelemento(e,edof,Rs,ex,ey,ez,xi,eta)

   Aux(1,:)=Rs(edof(e,2),:);
   Aux(2,:)=Rs(edof(e,3),:);
   Aux(3,:)=Rs(edof(e,4),:);
   Aux(4,:)=Rs(edof(e,5),:);
   A=Aux';
  

    N1=1/4*(xi-1).*(eta-1);
    N2=-1/4*(xi+1).*(eta-1);
    N3=1/4*(xi+1).*(eta+1); 
    N4=-1/4*(xi-1).*(eta+1);
    N=[N1 N2 N3 N4];
%--------------------------------------------------------------------
   %****Area do elemento****

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
Area=sqrt(detJxy.^2+detJyz.^2+detJzx.^2);
Area=sum(Area);