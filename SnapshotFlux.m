function [EdotSeries, Edot] = SnapshotFlux(a, e, p, iota_deg)

% calculate the geodesic
orbit = KerrGeodesic(a, e, p, iota_deg, 1e-10);

% calculate the cartesian metric (really a bunch of series representation
% of the cartesian metric's components expressed in spherical coordinates)
h = SeriesCartesianMetric(orbit);

% get hphx at a few places
mkn = [];
phi_deg = 0;
for theta_deg = [20 40 60]
    [hplus hcross] = Serieshphx(h, theta_deg, phi_deg);
    hcross.coefs = -i*hcross.coefs;
    hphx = SeriesAdd(hplus,hcross);
    hphx = TruncateSeries(hphx,1e-2);
    mkn = [mkn; hphx.mkn];
    mkn = unique(mkn,'rows');
end
[uniquemodes col] = size(mkn);

% compute the Teukolsky flux using executable temporary file
!rm TmpModeFile
for ell = 2:1:4
    parfor mkn_ind = 1:uniquemodes
        CallTeuk(a, p, e, iota_deg, ell, mkn(mkn_ind,:));
    end
end

% compute the energy flux series
EdotSeries = SeriesQuadEdot(h);

% the flux should be averaged over many wave cycles.  That average is 
% probably best done by setting the flux equal to the term in the series 
% with zero frequency.  Averaging analytically will fail because of this
% term, and averaging numerically approximates this term.
omega = EdotSeries.mkn * EdotSeries.freqs;
Edot = EdotSeries.coefs(omega == 0);

% strip off imaginary part, abort if it's not roundoff junk.
if abs(imag(Edot)/real(Edot)) > 1e-5
    error(['SnapshotFlux: Flux has significant imaginary part Edot = ' num2str(Edot) '.']);
end
Edot = real(Edot);
