% fiducial choice
S.a = 0.6;
S.e0 = 0.4;
S.p0 =1.035333114568974e+01;
S.iota0_deg = 2.998194158095480e+01;
S.r0 = S.p0 / ( 1 - S.e0);
S.theta0_deg = 90;
S.phi0_deg = -4.713797938634511e+06;
S.sign_rdot0 = 1;
S.sign_Tdot0 = -1;
S.M = 1.e6;
S.mu = 10.e0;
S.t0=-3.045358484207205e+06;
S.tspan = 3*2.5e7;
S.BigSteps = 16;
S.tol = 1e-6;
S.D = 1.e17;
S.thetasb_deg = 45;
S.phisb_deg = 150;
S.theta_k_deg = 60;
S.phi_k_deg = 200;
S.order = 'quadrupole';

% very simple 
% S.a = 0.1;
% S.e0 = 0.1;
% S.p0 =50;
% S.iota0_deg = 10;
% S.r0 = S.p0 / ( 1 - S.e0);
% S.theta0_deg = 90;
% S.phi0_deg = -4.713797938634511e+06;
% S.sign_rdot0 = 1;
% S.sign_Tdot0 = -1;
% S.M = 1.e6;
% S.mu = 10.e0;
% S.t0=0;
% S.tspan = 3e7;
% S.BigSteps = 16;
% S.tol = 1e-6;
% S.D = 1.e17;
% S.thetasb_deg = 45;
% S.phisb_deg = 150;
% S.theta_k_deg = 60;
% S.phi_k_deg = 200;
% S.order = 'quadrupole';

% the tail end of KITS: the EMRI that I also have Teukolsky spectra for 
% (the one in PRD, Drasco, 2009)
% S.a = 0.9;
% S.e0 = 6.984210e-01;
% S.p0 = 7.911813;
% S.iota0_deg = 4.292333e+01;
% S.r0 = S.p0 / (1 - S.e0);
% S.theta0_deg = 90;
% S.phi0_deg = 0;
% S.sign_rdot0 = -1;
% S.sign_Tdot0 = 1;
% S.M = 1.e6;
% S.mu = 10.e0;
% S.t0=9.0234e+06; % this is in M
% S.tspan = 0.99 * (9.439492e+07 - 4.44447567e+07); % this is in seconds
% S.BigSteps = 64;
% S.tol = 1e-5;
% S.D = 1.e17;
% S.thetasb_deg = 60;
% S.phisb_deg = 150;
% S.theta_k_deg = 60;
% S.phi_k_deg = 200;
% S.order = 'quadrupole';

% the "full length" KITS
% S.a = 0.9;
% S.e0=7.934560000000001e-01;
% S.p0=8.756240000000000e+00;
% S.iota0_deg=4.266550000000000e+01;
% S.r0 = S.p0 / (1 - S.e0);
% S.theta0_deg = 90;
% S.phi0_deg = 0;
% S.sign_rdot0 = -1;
% S.sign_Tdot0 = 1;
% S.M = 1.e6;
% S.mu = 10.e0;
% S.t0=0; % this is in M
% S.tspan = 0.99 * (9.439492e+07); % this is in seconds
% S.BigSteps = 32;
% S.tol = 1e-5;
% S.D = 1.e17;
% S.thetasb_deg = 60;
% S.phisb_deg = 150;
% S.theta_k_deg = 60;
% S.phi_k_deg = 200;
% S.order = 'quadrupole';

% % a Schwarzschild example similar to KITS
% S.a = 0.001;
% S.e0 = 6.984210e-01;
% S.p0 = 7.911813;
% S.iota0_deg = 4.292333e+01;
% S.r0 = S.p0 / (1 - S.e0);
% S.theta0_deg = 90;
% S.phi0_deg = 0;
% S.sign_rdot0 = -1;
% S.sign_Tdot0 = 1;
% S.M = 1.e6;
% S.mu = 10.e0;
% S.t0=9.0234e+06; % this is in M
% S.tspan = 0.07 * (9.439492e+07 - 4.44447567e+07); % this is in seconds
% S.BigSteps = 16;
% S.tol = 1e-5;
% S.D = 1.e17;
% S.thetasb_deg = 60;
% S.phisb_deg = 150;
% S.theta_k_deg = 60;
% S.phi_k_deg = 200;
% S.order = 'quadrupole';

% a circular example similar to KITS
% S.a = 0.9;
% S.e0 = 0.001;
% S.p0 = 7.911813;
% S.iota0_deg = 4.292333e+01;
% S.r0 = S.p0 / (1 - S.e0);
% S.theta0_deg = 90;
% S.phi0_deg = 0;
% S.sign_rdot0 = -1;
% S.sign_Tdot0 = 1;
% S.M = 1.e6;
% S.mu = 10.e0;
% S.t0=9.0234e+06; % this is in M
% S.tspan = 0.5 * (9.439492e+07 - 4.44447567e+07); % this is in seconds
% S.BigSteps = 16;
% S.tol = 1e-5;
% S.D = 1.e17;
% S.thetasb_deg = 60;
% S.phisb_deg = 150;
% S.theta_k_deg = 60;
% S.phi_k_deg = 200;
% S.order = 'quadrupole';



