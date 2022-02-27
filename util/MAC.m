function eta = MAC(a,b)
    a = a(:);
    b = b(:);
    eta = (a'*b)^2/((a'*a)*(b'*b));
end