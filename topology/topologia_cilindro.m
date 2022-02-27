function [coord,edof,conec,dof] = topologia_cilindro(raio,altura,Ne_arco,Ne_z)
    Ne = Ne_arco*Ne_z;
    Np = (Ne_arco+1)*(Ne_z+1);
    coord = zeros(Np,3);
    dteta = 2*pi/Ne_arco;
    p = 0;
    for z = 0:altura/Ne_z:altura
        for teta = 0:dteta:2*pi
            p = p+1;
            coord(p,1) = raio*sin(teta);
            coord(p,2) = raio*cos(teta);
            coord(p,3) = z;
        end
    end
    conec = zeros(Ne,4);
    for i = 1:Ne_z
        for j = 1:Ne_arco
            p1 = j + (Ne_arco+1)*(i-1);
            p2 = (j+1)+ (Ne_arco+1)*(i-1);
            p3 = (j+1)+ (Ne_arco+1)*i;
            p4 = j + (Ne_arco+1)*i;
            e = j + Ne_arco*(i-1);
            conec(e,:) = [p1 p2 p3 p4];
        end
    end       
    dof=(1:Np)';
    edof = zeros(Ne,5);
    edof(:,1) =(1:Ne)';
    edof(:,2:5) = conec;
 end