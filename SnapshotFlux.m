function [EdotSeries, Edot] = SnapshotFlux(a, e, p, iota_deg)

% calculate the geodesic
orbit = KerrGeodesic(a, e, p, iota_deg, 1e-10);

% calculate the cartesian metric (really a bunch of series representation
% of the cartesian metric's components expressed in spherical coordinates)
h = SeriesCartesianMetric(orbit);

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