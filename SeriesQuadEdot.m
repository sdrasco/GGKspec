function Edot=SeriesQuadEdot(h)
%
% Steve Drasco
%

% time derivatives of quadrupole metric (dh/dt)
dhxx = h.quad.dxx;
dhyy = h.quad.dyy;
dhzz = h.quad.dzz;
dhxy = h.quad.dxy;
dhxz = h.quad.dxz;
dhyz = h.quad.dyz;

% series for total flux integrated over sky, from mathematica calculation
%
% Edot = int_sky (hplus^2 + hcross^2) / (16 pi)
%
% Edot = dhxx^2/30 
%      + dhxy^2/10 
%      + dhxz^2/10 
%      - (dhxx dhyy)/30 
%      + dhyy^2/30 
%      + dhyz^2/10 
%      - (dhxx dhzz)/30 
%      - (dhyy dhzz)/30 
%      + dhzz^2/30
%
% This should also be averaged over time.

term = SeriesMultiply(dhxx,dhxx);
C = 1/30;
term.coefs = C*term.coefs;
Edot = term;

term = SeriesMultiply(dhxy,dhxy);
C = 1/10;
term.coefs = C*term.coefs;
Edot = SeriesAdd(Edot,term);

term = SeriesMultiply(dhxz,dhxz);
C = 1/10;
term.coefs = C*term.coefs;
Edot = SeriesAdd(Edot,term);

term = SeriesMultiply(dhxx,dhyy);
C = -1/30;
term.coefs = C*term.coefs;
Edot = SeriesAdd(Edot,term);

term = SeriesMultiply(dhyy,dhyy);
C = 1/30;
term.coefs = C*term.coefs;
Edot = SeriesAdd(Edot,term);

term = SeriesMultiply(dhyz,dhyz);
C = 1/10;
term.coefs = C*term.coefs;
Edot = SeriesAdd(Edot,term);

term = SeriesMultiply(dhxx,dhzz);
C = -1/30;
term.coefs = C*term.coefs;
Edot = SeriesAdd(Edot,term);

term = SeriesMultiply(dhyy,dhzz);
C = -1/30;
term.coefs = C*term.coefs;
Edot = SeriesAdd(Edot,term);

term = SeriesMultiply(dhzz,dhzz);
C = 1/30;
term.coefs = C*term.coefs;
Edot = SeriesAdd(Edot,term);








