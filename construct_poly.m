function P = construct_poly(ctrs,rbf)
%GOAL: This constructs the matrix P we typically need in all of our work.
%It is set to create polynomials of up to l=3, but depending on the order
%of the RBF, it restricts down to what is needed and leaves out extraneous
%values. 


%These polynomials come from the Wikipedia Table of Spherical Harmonics
%http://en.wikipedia.org/wiki/Table_of_spherical_harmonics
x = ctrs(:,1);
y = ctrs(:,2);
z = ctrs(:,3);

n = length(ctrs);
%Degree 0 Spherical Harmonic (constant)
P0 = .5*1/sqrt(pi)*ones(n,1);

%Degree 1 Spherical Harmonic 
P1 =sqrt(3/(4*pi))*[x y z];

%Degree 2 Spherical Harmonics
P2 = .25*sqrt(5/pi)*[-x.^2-y.^2+2*z.^2  2*sqrt(3)*y.*z  2*sqrt(3)*z.*x   2*sqrt(3)*x.*y   sqrt(3)*(x.^2-y.^2)];

%Degree 3 Spherical Harmonics
P3 = [z.*(2*z.^2-3*x.^2-3*y.^2)   (3*x.^2-y.^2).*y   ...
    (x.^2-3*y.^2).*x   (x.^2-y.^2).*z  x.*y.*z ...
    y.*(4*z.^2-x.^2-y.^2)   x.*(4*z.^2-x.^2-y.^2)];

P = [P0 P1 P2 P3];
%Only choose necessary polynomials
P = P(:,1:(rbf+1)^2);

