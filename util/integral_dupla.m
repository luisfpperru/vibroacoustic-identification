function I = integral_dupla(hx,hy,F)
    nc = size(F,2);
    Ix = zeros(nc,1);
    for j = 1:nc
        Ix(j) = Simpson3(hx,F(:,j)); 
    end
    I = Simpson3(hy,Ix); 
end