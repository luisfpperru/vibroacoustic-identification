% Integral pela regra de Simpson de 3 pontos

function I = Simpson3(h,fx)
    fx = fx(:);
    n = size(fx,1);
    soma = 0;
    for i = 2:n-1
        if (-1)^i == 1
                soma = soma + 4*fx(i);
        end
        if (-1)^i == -1
                 soma = soma + 2*fx(i);
        end
    end
    I = h/3*(fx(1)+fx(n)+soma);
end