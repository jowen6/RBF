function [A,rc] =  covariant_stiffness(alpha,beta,ctrs,ep,qwe)

[Nchi1,Nchi2,Nchi3] = phi_grad(ctrs,ep);

NC1 = Nchi1*alpha;
NC2 = Nchi2*alpha;
NC3 = Nchi3*alpha;
[Px, Py, Pz] = poly_grad(ep);
covchi1 = NC1 + Px*beta;
covchi2 = NC2 + Py*beta;
covchi3 = NC3 + Pz*beta;

A = zeros(length(ctrs),length(ctrs));
for xi=1:length(ctrs)
    xi;
    for eta=xi:length(ctrs)
        A(xi,eta) = dot(covchi1(:,xi).*covchi1(:,eta)+covchi2(:,xi).*covchi2(:,eta)+covchi3(:,xi).*covchi3(:,eta),qwe);
        A(eta,xi) = A(xi,eta);
    end
end
rc = rcond(A);