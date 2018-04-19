function [Nphi1,Nphi2,Nphi3] = phi_grad(ctrs,ep)

n = length(ep);
x = ep(:,1);
y = ep(:,2);
z = ep(:,3);

theta = acos(z);

phi = atan2(y,x);

theta_hat =[ cos(theta).*cos(phi) cos(theta).*sin(phi) -sin(theta)];

phi_hat = [-sin(phi) cos(phi) zeros(n,1)];


gradx = repmat(cos(theta).*cos(phi),1,3).*theta_hat -repmat(sin(phi),1,3).*phi_hat;

grady = repmat(cos(theta).*sin(phi),1,3).*theta_hat +repmat(cos(phi),1,3).*phi_hat;

gradz = repmat(-sin(theta),1,3).*theta_hat;

phi2e = -2*(1-ep*ctrs').*log(1-ep*ctrs')+ep*ctrs'-1;
phi2e(isnan(phi2e))=0;

comp_1 = ctrs(:,1)*gradx(:,1)' + ctrs(:,2)*grady(:,1)' + ctrs(:,3)*gradz(:,1)';

comp_2 = ctrs(:,1)*gradx(:,2)' + ctrs(:,2)*grady(:,2)' + ctrs(:,3)*gradz(:,2)';

comp_3 = ctrs(:,1)*gradx(:,3)' + ctrs(:,2)*grady(:,3)' + ctrs(:,3)*gradz(:,3)';

Nphi1 = phi2e.*comp_1';
Nphi2 = phi2e.*comp_2';
Nphi3 = phi2e.*comp_3';