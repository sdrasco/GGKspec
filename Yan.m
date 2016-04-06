
% fiducial choice
S.a = 0.7;
S.e0 = 0.86;
S.p0 = 9;
S.iota0_deg = 50;
S.r0 = S.p0 / ( 1 - S.e0);
S.theta0_deg = 90;
S.phi0_deg = 0;
S.sign_rdot0 = 1;
S.sign_Tdot0 = -1;
S.M = 1e6;
S.mu = 1.4;
S.t0=0;
S.tspan = 3e7;
S.BigSteps = 10;
S.tol = 1e-3;
S.D = 1.e17;
S.thetasb_deg = 45;
S.phisb_deg = 150;
S.theta_k_deg = 60;
S.phi_k_deg = 200;

% get the initial waveform data structure
waveform = InitializeWaveformSeries(S);