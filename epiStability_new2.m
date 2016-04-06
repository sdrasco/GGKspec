%8/26/10 This is for testing Steve's new version of code
%
%Jan., 2010 :In this version we try 4-pt diff formulae as well as 2-pt ones
%
%in this version both A and B version of derivs are splined
%directly to same t vector
delta = 1.e-9;
%eps = 1.e-5;
eps = 1.e-8;


%first take finite derivative at pt A
muA = 300.e0;
%Curt choosing loarge muA to compensate for shorter integration time:
%muA =1.e2;
% load sample parameters
% the "full length" KITS
S.a = 0.9;
%S.e0=7.934560000000001e-01;
%Curt lowered e0
S.e0=3.934560000000001e-01;
S.p0=8.756240000000000e+00;
S.iota0_deg=4.266550000000000e+01;
S.r0 = S.p0 / (1 - S.e0);
S.theta0_deg = 90;
S.phi0_deg = 0;
S.sign_rdot0 = -1;
S.sign_Tdot0 = 1;
S.M = 1.e6;
S.mu = muA;
S.t0=0; % this is in M
%S.tspan = 0.99 * (9.439492e+07); % this is in seconds
S.tspan = 0.99 * (9.439492e+05); % Curt divided tspan by 10 to save time in testing
S.BigSteps = 16;
%S.BigSteps = 4; %Curt changed this to 4 (from 32) to save time in testing
S.tol = 1e-10;
S.D = 1.e17;
S.thetasb_deg = 60;
S.phisb_deg = 150;
S.theta_k_deg = 60;
S.phi_k_deg = 200;
S.order = 'quadrupole';

% steve's edits to make parameters more like those from first generation of
% stability tests
S.a = 0.6;
S.e0 = 3.177102857980181e-01;
S.p0 =1.035333114568974e+01;
S.iota0_deg = 2.998194158095480e+01;
S.M = 1.e6;
S.mu = 10.e0;
S.r0 = S.p0 / (1 - S.e0);
S.theta0_deg = 90;
S.tspan = 3e6;

%now take finite deriv wrt mu
mup = muA*(1.0 + eps);
mum = muA*(1.0 - eps);
TrueEps = 0.5*(mup - mum)/muA;

S.mu = mup
waveformpA =  InitializeWaveformSeries(S);
wspA = SplineWaveformSeries(waveformpA);
[hpA wpA tp] = StitchMode(wspA,2,0,0,10^4);
tpA = tp;

S.mu = mum
waveformmA =  InitializeWaveformSeries(S);
wsmA = SplineWaveformSeries(waveformmA);
[hmA wmA tm] = StitchMode(wsmA,2,0,0,10^4);
tmA = tm;

dhA = (hpA - hmA)/(2.0*TrueEps);
dwA = (wpA - wmA)/(2.0*TrueEps);
display('finished pt A');

clear mup ; clear mum; clear TrueEps;
%Curt: rem to also clear hp,wp after some checks
%now repeat at pt. B
muB = muA*(1.0 + delta);

%now take finite deriv wrt mu
mup = muB*(1.0 + eps);
mum = muB*(1.0 - eps);
TrueEps = 0.5*(mup - mum)/muB;

S.mu = mup
waveformpB =  InitializeWaveformSeries(S);
wspB = SplineWaveformSeries(waveformpB);
[hpB wpB tp] = StitchMode(wspB,2,0,0,10^4);
tpB = tp;

S.mu = mum
waveformmB =  InitializeWaveformSeries(S);
wsmB = SplineWaveformSeries(waveformmB);
[hmB wmB tm] = StitchMode(wsmB,2,0,0,10^4);
tmB = tm;

dhB = (hpB - hmB)/(2.0*TrueEps);
dwB = (wpB - wmB)/(2.0*TrueEps);


display('finished pt B');
figure;
plot(tp,dwB-dwA);
set(gca,'fontsize',16);
title('dwB-dwA');
figure;
plot(tp,real(dhB-dhA)./real(dhA));
set(gca,'fontsize',16);
title('real(dhB-dhA)./real(dhA)');
figure;
plot(tp,imag(dhB-dhA)/max(real(dhA)) );
set(gca,'fontsize',16);
title('imag(dhB-dhA)/max(real(dhA))');

%
semilogy(tp,abs(real(dhB-dhA))./abs(real(dhA)));set(gca,'fontsize',16);
semilogy(tp,abs(dwB-dwA)./abs(dwA));set(gca,'fontsize',16);
%
% Steve's block for looking at derivatives of other parameters with respect
% to mu.  As is, it's the derivative of E with respect to mu.  Change the
% four places where E appears to see another derivative.
%
close all;
Xtstart = 1.01*max([wspA.t(2) wspB.t(2)]);
Xtend = 0.99*min([wspA.t(end) wspB.t(end)]);
Xt = linspace(Xtstart, Xtend,1e4);
XpA = ppval(wspA.E_pp,Xt);
XmA = ppval(wsmA.E_pp,Xt);
dXA = (XpA - XmA)/(2.0*TrueEps);
XpB = ppval(wspB.E_pp,Xt);
XmB = ppval(wsmB.E_pp,Xt);
dXB = (XpB - XmB)/(2.0*TrueEps);
semilogy(Xt,abs(dXA-dXB)./abs(dXA));set(gca,'fontsize',16);


display('stopped');
%Curt: rem there are lots of examples of how to plot this stuff old
%directory in epiStability3.m
%%%
