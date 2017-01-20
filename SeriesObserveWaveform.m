function [hplus hcross]=SeriesObserveWaveform(hin, r, theta, phi, order)
%
% [hplus hcross]=SeriesObserveWaveform(waveform, r, theta, phi)
%
% Series version of ObserveWaveform (quadrupole only for now)
%
%      r: radial distance in units of M
%  theta: polar angle in degrees
%    phi: azimuthal angle in degrees
%
% See also OBSERVEWAVEFORM
%
% By: Steve Drasco
%

% convert to radians
theta = theta * pi/180;
phi = phi * pi/180;

if strcmp(order,'quadrupole')

    % use only the quadrupole formula
    h.xx = hin.quad.xx;
    h.yy = hin.quad.yy;
    h.zz = hin.quad.zz;
    h.xy = hin.quad.xy;
    h.xz = hin.quad.xz;
    h.yz = hin.quad.yz;

elseif strcmp(order,'octupole') % WARNING: doesn't pass tests

    % compute unit vector from source to observer;
    nx = 2*sin(theta) * cos(phi);
    ny = 2*sin(theta) * sin(phi);
    nz = 2*cos(theta);

    % project the octupole term, and sum with quadrupole term
    % xx
    h.xx = hin.quad.xx;
    nxterm = hin.oct.xxx;
    nxterm.coefs = nxterm.coefs * nx;
    nyterm = hin.oct.yxx;
    nyterm.coefs = nyterm.coefs * ny;
    nzterm = hin.oct.zxx;
    nzterm.coefs = nzterm.coefs * nz;
    h.xx = SeriesAdd(h.xx,nxterm);
    h.xx = SeriesAdd(h.xx,nyterm);
    h.xx = SeriesAdd(h.xx,nzterm);
    % yy
    h.yy = hin.quad.yy;
    nxterm = hin.oct.xyy;
    nxterm.coefs = nxterm.coefs * nx;
    nyterm = hin.oct.yyy;
    nyterm.coefs = nyterm.coefs * ny;
    nzterm = hin.oct.zyy;
    nzterm.coefs = nzterm.coefs * nz;
    h.yy = SeriesAdd(h.yy,nxterm);
    h.yy = SeriesAdd(h.yy,nyterm);
    h.yy = SeriesAdd(h.yy,nzterm);
    % zz
    h.zz = hin.quad.zz;
    nxterm = hin.oct.xzz;
    nxterm.coefs = nxterm.coefs * nx;
    nyterm = hin.oct.yzz;
    nyterm.coefs = nyterm.coefs * ny;
    nzterm = hin.oct.zzz;
    nzterm.coefs = nzterm.coefs * nz;
    h.zz = SeriesAdd(h.zz,nxterm);
    h.zz = SeriesAdd(h.zz,nyterm);
    h.zz = SeriesAdd(h.zz,nzterm);
    % xy
    h.xy = hin.quad.xy;
    nxterm = hin.oct.xxy;
    nxterm.coefs = nxterm.coefs * nx;
    nyterm = hin.oct.yxy;
    nyterm.coefs = nyterm.coefs * ny;
    nzterm = hin.oct.zxy;
    nzterm.coefs = nzterm.coefs * nz;
    h.xy = SeriesAdd(h.xy,nxterm);
    h.xy = SeriesAdd(h.xy,nyterm);
    h.xy = SeriesAdd(h.xy,nzterm);
    % xz
    h.xz = hin.quad.xz;
    nxterm = hin.oct.xxz;
    nxterm.coefs = nxterm.coefs * nx;
    nyterm = hin.oct.yxz;
    nyterm.coefs = nyterm.coefs * ny;
    nzterm = hin.oct.zxz;
    nzterm.coefs = nzterm.coefs * nz;
    h.xz = SeriesAdd(h.xz,nxterm);
    h.xz = SeriesAdd(h.xz,nyterm);
    h.xz = SeriesAdd(h.xz,nzterm);
    % yz 
    h.yz = hin.quad.yz;
    nxterm = hin.oct.xyz;
    nxterm.coefs = nxterm.coefs * nx;
    nyterm = hin.oct.yyz;
    nyterm.coefs = nyterm.coefs * ny;
    nzterm = hin.oct.zyz;
    nzterm.coefs = nzterm.coefs * nz;
    h.yz = SeriesAdd(h.yz,nxterm);
    h.yz = SeriesAdd(h.yz,nyterm);
    h.yz = SeriesAdd(h.yz,nzterm);

end

% we'll use these a bunch of times, and they are a pain to retype each time
cts = cos(theta)^2;
sts = sin(theta)^2;
cps = cos(phi)^2;
sps = sin(phi)^2;

% theta-theta component
%hTT = cts * (h.xx*cps + h.xy*sin(2*phi) + h.yy*sps ) ...
%    + h.zz*sts - sin(2*theta)*(h.xz*cos(phi) + h.yz *sin(phi));
term1 = h.xx;
term1.coefs = term1.coefs * cps * cts;
term2 = h.xy;
term2.coefs = term2.coefs * sin(2*phi) * cts ;
term3 = h.yy;
term3.coefs = term3.coefs * sps * cts;
term4 = h.zz;
term4.coefs = term4.coefs * sts;
term5 = h.xz;
term5.coefs = -term5.coefs * cos(phi) * sin(2*theta);
term6 = h.yz;
term6.coefs = -term6.coefs * sin(phi) * sin(2*theta);
hTT = SeriesAdd(term1,term2);
hTT = SeriesAdd(hTT,term3);
hTT = SeriesAdd(hTT,term4);
hTT = SeriesAdd(hTT,term5);
hTT = SeriesAdd(hTT,term6);

% theta-phi component
%hTP = cos(theta)*(h.xy*cos(2*phi) + 0.5*h.yy*sin(2*phi) - 0.5*h.xx*sin(2*phi)) ...
%    + sin(theta)*(h.xz*sin(phi) - h.yz*cos(phi));
term1 = h.xy;
term1.coefs = term1.coefs * cos(2*phi) * cos(theta);
term2 = h.yy;
term2.coefs = term2.coefs * 0.5 * sin(2*phi) * cos(theta);
term3 = h.xx;
term3.coefs = -term3.coefs * 0.5 * sin(2*phi) * cos(theta);
term4 = h.xz;
term4.coefs = term4.coefs * sin(phi) * sin(theta);
term5 = h.yz;
term5.coefs = -term5.coefs * cos(phi) * sin(theta);
hTP = SeriesAdd(term1,term2);
hTP = SeriesAdd(hTP,term3);
hTP = SeriesAdd(hTP,term4);
hTP = SeriesAdd(hTP,term5);

% phi-phi component
%hPP = h.xx*sps - h.xy*sin(2*phi) + h.yy*cps;
term1 = h.xx;
term1.coefs = term1.coefs * sps;
term2 = h.xy;
term2.coefs = -term2.coefs * sin(2*phi);
term3 = h.yy;
term3.coefs = term3.coefs * cps;
hPP = SeriesAdd(term1,term2);
hPP = SeriesAdd(hPP,term3);

% plus and cross
hcross = hTP;
hcross.coefs = hcross.coefs / r;
term1 = hTT;
term1.coefs = term1.coefs * 0.5;
term2 = hPP;
term2.coefs = -0.5 * term2.coefs;
hplus = SeriesAdd(term1,term2);
hplus.coefs = hplus.coefs / r;