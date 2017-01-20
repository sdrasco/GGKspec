function [hplus hcross] = Serieshphx(h,theta_deg,phi_deg)
%
% Steve Drasco
%
% Compute hplus and hcross from cartesian metric perturbation (h) using Eqs.
% (11.46) in Poisson and Will.

% convert angles to radians
theta = theta_deg * pi/180;
phi = phi_deg * pi/180;

% h_plus
%
% xx term
S = h.quad.xx;
C = -0.25*sin(theta)*sin(theta) + 0.25*(1+cos(theta)*cos(theta))*cos(2*phi);
S.coefs = C*S.coefs;
hplus = S;
% yy term
S = h.quad.yy;
C = -0.25*sin(theta)*sin(theta) - 0.25*(1+cos(theta)*cos(theta))*cos(2*phi);
S.coefs = C*S.coefs;
hplus = SeriesAdd(hplus,S);
% xyterm
S = h.quad.xy;
C = 0.5*(1+cos(theta)*cos(theta))*sin(2*phi);
S.coefs = C*S.coefs;
hplus = SeriesAdd(hplus,S);
% xzterm
S = h.quad.xz;
C = -sin(theta)*cos(theta)*cos(phi);
S.coefs = C*S.coefs;
hplus = SeriesAdd(hplus,S);
% yzterm
S = h.quad.yz;
C = -sin(theta)*cos(theta)*sin(phi);
S.coefs = C*S.coefs;
hplus = SeriesAdd(hplus,S);
% zzterm
S = h.quad.zz;
C = -sin(theta)*cos(phi);
S.coefs = C*S.coefs;
hplus = SeriesAdd(hplus,S);

% h_cross
% 
% xx term
S = h.quad.xx;
C = -0.5*cos(theta)*sin(2*phi);
S.coefs = C*S.coefs;
hcross = S;
% yy term
S = h.quad.yy;
C = 0.5*cos(theta)*sin(2*phi);
S.coefs = C*S.coefs;
hcross = SeriesAdd(hcross,S);
% xy term
S = h.quad.xy;
C = cos(theta)*cos(2*phi);
S.coefs = C*S.coefs;
hcross = SeriesAdd(hcross,S);
% xz term
S = h.quad.xz;
C = sin(theta)*sin(phi);
S.coefs = C*S.coefs;
hcross = SeriesAdd(hcross,S);
% yz term
S = h.quad.yz;
C = -sin(theta)*cos(phi);
S.coefs = C*S.coefs;
hcross = SeriesAdd(hcross,S);


