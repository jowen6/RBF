function [PX, PY, PZ] = poly_grad(ctrs)

x = ctrs(:,1);
y = ctrs(:,2);
z = ctrs(:,3);

n = length(ctrs);




theta = acos(z);

phi = atan2(y,x);

theta_hat =[ cos(theta).*cos(phi) cos(theta).*sin(phi) -sin(theta)];

phi_hat = [-sin(phi) cos(phi) zeros(n,1)];


gradx = repmat(cos(theta).*cos(phi),1,3).*theta_hat -repmat(sin(phi),1,3).*phi_hat;

grady = repmat(cos(theta).*sin(phi),1,3).*theta_hat +repmat(cos(phi),1,3).*phi_hat;

gradz = repmat(-sin(theta),1,3).*theta_hat;

rx = [x x x];
ry = [y y y];
rz = [z z z];
P02 =  .25*sqrt(5/pi)*(-2*rx.*gradx - 2*ry.*grady +4*rz.*gradz);

Pyz = .5*sqrt(15/pi)*(ry.*gradz + rz.*grady);
Pzx = .5*sqrt(15/pi)*(rz.*gradx + rx.*gradz);
Pxy = .5*sqrt(15/pi)*(rx.*grady + ry.*gradx);

Px2y2 = .25*sqrt(15/pi)*(2*rx.*gradx - 2*ry.*grady);
c1 =  sqrt(3/(4*pi));

PX =[zeros(n,1) c1*gradx(:,1) c1*grady(:,1) c1*gradz(:,1) ...
    P02(:,1) Pyz(:,1) Pzx(:,1) Pxy(:,1) Px2y2(:,1)];

PY =[zeros(n,1) c1*gradx(:,2) c1*grady(:,2) c1*gradz(:,2) ...
    P02(:,2) Pyz(:,2) Pzx(:,2) Pxy(:,2) Px2y2(:,2)];

PZ =[zeros(n,1) c1*gradx(:,3) c1*grady(:,3) c1*gradz(:,3) ...
    P02(:,3) Pyz(:,3) Pzx(:,3) Pxy(:,3) Px2y2(:,3)];