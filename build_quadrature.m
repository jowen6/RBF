function [ep,qwe] = build_quadrature(j)

epstr = sprintf('%04d',(j+1)^2);
epstr2 = ['me' num2str(j) '.' epstr];
ep = load(epstr2);
ep = ep(:,1:3);
disp('Building Quadrature Weights')
ne = length(ep);
phiep = (1-ep*ep').*log(1-ep*ep');
phiep(isnan(phiep))=0;
phiep = real(phiep);
PE = construct_poly(ep,1);
p = min(size(PE));
ecoefs = [phiep PE; PE' zeros(p,p)]\[eye(ne,ne);zeros(p,ne)];
clear phiep


%Quadrature weights proportional to beta_0 coefficients
qwe = 2*sqrt(pi)*ecoefs(ne+1,:)';
clear ecoefs