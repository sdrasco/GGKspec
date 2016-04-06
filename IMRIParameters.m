% fiducial choice
S.a = 0.5;
S.e0 = 0.8;
S.p0 =10;
S.iota0_deg = 33;
S.r0 = S.p0 / ( 1 - S.e0);
S.theta0_deg = 90;
S.phi0_deg = 0;
S.sign_rdot0 = 1;
S.sign_Tdot0 = -1;
S.M = 250;
S.mu = 1.4;
S.t0=0;
S.tspan = 35.15;
S.BigSteps = 8;
S.tol = 1e-3;

% these last parameters tell where where the binary is in the sky, and how
% far it is (S.D)
S.D = 1.e17;
S.thetasb_deg = 45;
S.phisb_deg = 150;
S.theta_k_deg = 60;
S.phi_k_deg = 200;
S.order = 'quadrupole';






