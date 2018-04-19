function [ctrs, coefs] = build_ctrs(j)

    numstr2 = sprintf('%04d',(j+1)^2);
    ctr_str = ['me' num2str(j) '.' numstr2];
    ctrs = load(ctr_str);
    ctrs = ctrs(:,1:3);
    n = length(ctrs);
    phi = (1-ctrs*ctrs').^2.*log(1-ctrs*ctrs');
    phi(1:n+1:n^2)=0;
    P = construct_poly(ctrs,2);
    p = min(size(P));
    coefs = [phi P; P' zeros(p,p)]\[eye(n,n);zeros(p,n)];