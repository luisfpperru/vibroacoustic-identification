function eta = Spearman(x,y)
    x = x(:);
    y = y(:);
    [~,orderx] = sort(x);
    [~,ordery] = sort(y);
    d = orderx - ordery;
    n = size (d, 1);
    eta = 1 - 6*(d'*d)/(n*(n^2 -1));
end