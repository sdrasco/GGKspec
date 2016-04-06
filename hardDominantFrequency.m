%function Omega_in_Hz_pp = hardDominantFrequency()
function [Omega_in_Hz mkn Omega_in_Hz_pp] = hardDominantFrequency()

% fiducial choice
S.a = 0.9;
S.e0 = 0.5;
S.p0 =6;
S.iota0_deg = 50;
S.r0 = S.p0 / ( 1 - S.e0);
S.theta0_deg = 90;
S.phi0_deg = 0;
S.sign_rdot0 = 1;
S.sign_Tdot0 = -1;
S.M = 3000;
S.mu = 1.4;
S.t0=0;
S.tspan = 643.0725;
S.BigSteps = 64;
S.tol = 1e-3;
S.D = 1.e17;
S.thetasb_deg = 45;
S.phisb_deg = 150;
S.theta_k_deg = 60;
S.phi_k_deg = 200;

% get the initial waveform data structure
%waveform = InitializeWaveformSeries(S);
%save hardIMRIwaveform.mat waveform
load hardIMRIwaveform.mat

% find dominant frequency at each BigStep
DominantFrequency = zeros(S.BigSteps,1);
mkn = zeros(S.BigSteps,3);
for ii=1:S.BigSteps
    p = abs(waveform.H{ii}.coefs);
    ind = find(p == max(p));
    mkn(ii,:) = waveform.H{ii}.mkn(ind,:);
    DominantFrequency(ii) = mkn(ii,:) * waveform.H{ii}.freqs;
end

% switch to physicaly meaningful units, and return a spline of omega(t)
SecPerMsun = 4.9255e-6;
SecPerM = waveform.M * SecPerMsun;
Omega_in_Hz = DominantFrequency / SecPerM;
Omega_in_Hz_pp = spline(waveform.t,Omega_in_Hz);
