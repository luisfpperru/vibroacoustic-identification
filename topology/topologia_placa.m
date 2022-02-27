function [coord,edof,conec,dof] = topologia_placa(Lx,Ly,Ne_x,Ne_y)
    Ne = Ne_x*Ne_y;
    Np = (Ne_x+1)*(Ne_y+1);
    coord = zeros(Np,3);
    p = 0;
    for y = 0:Ly/Ne_y:Ly
        for x = 0:Lx/Ne_x:Lx
            p = p+1;
            coord(p,1) = x;
            coord(p,2) = y;
            coord(p,3) = 0;
        end
    end
    conec = zeros(Ne,4);
    for i = 1:Ne_y
        for j = 1:Ne_x
            p1 = j + (Ne_x+1)*(i-1);
            p2 = (j+1)+ (Ne_x+1)*(i-1);
            p3 = (j+1)+ (Ne_x+1)*i;
            p4 = j + (Ne_x+1)*i;
            e = j + Ne_x*(i-1);
            conec(e,:) = [p1 p2 p3 p4];
        end
    end       
    dof=(1:Np)';
    edof = zeros(Ne,5);
    edof(:,1) =(1:Ne)';
    edof(:,2:5) = conec;
 end