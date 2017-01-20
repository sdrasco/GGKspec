function waveform = InitializeWaveformSeries(S)
% 
% waveform = InitializeWaveform_epiELQ(S);
%
% A series version of InitializeWaveform_epiELQ.m.
%
% Input structure has the following fields:
%
%          a: dimensionless spin of larger black hole [0, 1]
%         e0: initial eccentricity (0, 1)
%         p0: initial semilatus rectum [0 inf]
%  iota0_deg: initial inclination in degrees [0 180] 
%             defined as iota = 90 - sgn(L) theta_min
%         t0: t(lambda=0) in M             [-inf inf]
%         r0: r(lambda=0) in M             [p/(1+e)  p/(1-e)]
% theta0_deg: theta(lambda=0)  in degrees [theta_min  pi-theta_min], where
%             theta_min = 90 - iota0_deg for prograde orbits (iota < 90),
%             theta_min = iota0_deg - 90 for retrograde orbits (iota > 90).
%   phi0_deg: phi(lambda=0) in degrees    [0  360] 
% sign_rdot0: sign of initial radial velocity [-1 or 1]
% sign_Tdot0: sign of initial theta velocity [-1 or 1]
%          M: mass of the large black hole, in solar masses (mu inf]
%         mu: mass of the small black hole, in solar masses [0 M)
%      tspan: final time minus initial time, in seconds (-inf inf)
% SmallSteps: number of time steps for waveform [1 inf]
%   BigSteps: number of time steps for the trajectory of the orbital
%             elements [3 inf]
%        tol: accuracy to which we will compute the orbital paramaters of 
%             the snapshots (including the coefficients in their series 
%             expansions).
%     coords: a string telling whether to use spherical (coords =
%             'spherical') or spheroidal (coords = 'spheroidal') 
%             coordinates.
%
% See also INITIALIZEWAVEFORM_EPIELQ PEIT
% 
% By: Steve Drasco
% 

% strip input parameters from structure
% 
% This block also functions as a check that input structure has all the
% necessary fields
a=S.a;
e0=S.e0;
p0=S.p0;
iota0_deg=S.iota0_deg;
t0=S.t0;
r0=S.r0;
theta0_deg=S.theta0_deg;
phi0_deg=S.phi0_deg;
sign_rdot0=S.sign_rdot0;
sign_Tdot0=S.sign_Tdot0;
M=S.M;
mu=S.mu;
tspan=S.tspan;
BigSteps=S.BigSteps;
tol=S.tol;

% start clock to keep track of computational cost
InitialCPUTime = cputime;

% make sure that inputs are reasonable
if e0 < 0 || e0 >= 1
    error('We require 0 < e0 < 1.');
end
if p0 < 0
    error('We require a positive p0.');
end
if iota0_deg > 180 || iota0_deg < 0
    error('We require 0 < iota0_deg < 180.');
end
if M < 0
    error('We require M > 0.');
end
if mu < 0 || mu > M
    error('We require 0 < mu < M.');
end
if BigSteps < 3
    error('We require 2 < BigSteps');
end
if abs(sign_rdot0) ~= 1
    error('sign_rdot0 must be 1 or -1.');
end
if abs(sign_Tdot0) ~= 1
    error('sign_Tdot0 must be 1 or -1.');
end
if r0 < p0/(1 + e0) || r0 > p0/(1 - e0)
    error('InitializeWaveform: r0 is out of bounds.');
end
if iota0_deg < 90
  theta_min = 90 - iota0_deg;
else
  theta_min = iota0_deg - 90;
end
if theta0_deg < theta_min || theta0_deg > (180-theta_min)
    error('InitializeWaveform: theta0_deg is out of bounds.');
end

% put the input parameters into the output structure
waveform.a = a;
waveform.e0 = e0;
waveform.p0 = p0;
waveform.iota0_deg = iota0_deg;
waveform.iota_hughes_deg0 = HughesIota(iota0_deg, a, e0, p0);
waveform.t0 = t0;
waveform.r0 = r0;
waveform.theta0_deg = theta0_deg;
waveform.phi0_deg = phi0_deg;
waveform.sign_rdot0=sign_rdot0;
waveform.sign_Tdot0=sign_Tdot0;
waveform.M = M;
waveform.mu = mu;
waveform.BigSteps = BigSteps;
waveform.tol = tol;
waveform.SecPerMsun = 4.9255e-6;
waveform.SecPerM = waveform.SecPerMsun * waveform.M;
waveform.SecPermu = waveform.SecPerMsun * waveform.mu;
[waveform.thetaE_deg waveform.phiE_deg] = EarthSkyPosition(S);
if waveform.phiE_deg ~= 0
    error('Your phi_E is not zero.  We are currently using Curt''s formula for the waveform, which requires phi_E = 0.');
end

% evolve the principle orbital elements: 
% make a trajectory e(t) p(t) iota(t), using a big time step
[t, p, e, iota_hughes_rad, E, L, Q] = peitELQ(p0,e0,waveform.iota_hughes_deg0*pi/180,...
    t0*waveform.SecPerM,t0*waveform.SecPerM + tspan,BigSteps,M,a*M,mu);
waveform.t = t;
waveform.p = p;
waveform.e = e;
waveform.E = E / waveform.mu;
waveform.L = L / (waveform.mu * waveform.M);
waveform.Q = Q / (waveform.mu * waveform.M)^2;

% translate to geometric definition of inclination
iota_deg = zeros(size(iota_hughes_rad));
for i=1:length(iota_hughes_rad)
    iota_deg(i)=GeometricIota(iota_hughes_rad(i)*180/pi, a, e(i), p(i), 1e-11);
end
waveform.iota_deg = iota_deg;

% translate t into units of M
waveform.t_in_M = t / waveform.SecPerM; 
waveform.t_in_years = t / 31556926;

% for each time, we generate a snapshot
for j = 1:waveform.BigSteps
    display(['working on geodesic number ' num2str(j) ' / ' num2str(waveform.BigSteps)]);
    waveform.orbit{j} = KerrGeodesic(a, e(j), p(j), iota_deg(j), tol, waveform.E(j), waveform.L(j), waveform.Q(j));
    snapshot = SnapshotSeries(waveform.orbit{j});
    waveform.H{j} = ObserveSnapshot(snapshot,waveform.thetaE_deg);
    waveform.H{j}.coefs = waveform.H{j}.coefs * (waveform.mu * waveform.SecPerMsun) / S.D;
end


% log the computational cost of this job
waveform.CPUsec = cputime - InitialCPUTime;

