function snapshot = SnapshotSeries(orbit)
%
%
% Steve Drasco
%

% strip tol from orbit
tol = orbit.tol;

% we need r^2, sin(theta), and cos(theta) many times
r2 = SeriesMultiply(orbit.r,orbit.r);
cosT = SeriesRe(orbit.eiT);
sinT = SeriesIm(orbit.eiT);

% first term
r2cosT = SeriesMultiply(r2,cosT);
r2cos2T = SeriesMultiply(r2cosT,cosT);
A = r2cos2T;
A.coefs = 1.5 * A.coefs;
B = r2;
B.coefs = -0.5 * B.coefs;
term1 = SeriesAdd(A,B);
term1 = SeriesDifferentiate(term1);
term1 = SeriesDifferentiate(term1);

% second term
r2cosTsinT = SeriesMultiply(r2cosT,sinT);
term2 = SeriesMultiply(r2cosTsinT,orbit.eiphi);
term2 = SeriesDifferentiate(term2);
term2 = SeriesDifferentiate(term2);

% third term
term3 = SeriesConj(term2);

% fourth term
e2iphi = SeriesMultiply(orbit.eiphi,orbit.eiphi);
%cos2phi = SeriesRe(e2iphi);
sin2T = SeriesMultiply(sinT,sinT);
r2sin2T = SeriesMultiply(r2,sin2T);
r2sin2Te2iphi = SeriesMultiply(r2sin2T,e2iphi);
term4 = SeriesRe(r2sin2Te2iphi);
term4 = SeriesDifferentiate(term4);
term4 = SeriesDifferentiate(term4);

% fifth term
term5 = SeriesIm(r2sin2Te2iphi);
term5.coefs = (1i) * term5.coefs;
term5 = SeriesDifferentiate(term5);
term5 = SeriesDifferentiate(term5);

% collect results
snapshot.term1 = term1;
snapshot.term2 = term2;
snapshot.term3 = term3;
snapshot.term4 = term4;
snapshot.term5 = term5;

