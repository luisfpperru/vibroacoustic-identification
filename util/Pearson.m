function eta = Pearson(x,y)
    eta = (mean(x.*y) - mean(x)*mean(y))/(sqrt(mean(x.^2) - mean(x)^2)*sqrt(mean(y.^2) - mean(y)^2));
end