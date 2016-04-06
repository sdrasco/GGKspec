function waveform = KITSconstructor(theta_rad)
% 
% This is a special routine which uses the Teukoslky data set KITS to
% make a structure which has fields that are similar to the fields in the
% structures produced by InitializeWaveformSeries.  The differences are:
%
%  H: has the same format, but is made from the Teukolsky snapshots
% 
% By: Steve Drasco
% 

% start clock to keep track of computational cost
InitialCPUTime = cputime;

% load the SnapshotMatrix, which gives the trajectory of orbital parameters
% used for KITS
SM = load('/data/KITS/survey/SnapshotMatrix');

% read in the simple orbital parameters
waveform.t = SM(:,1);
waveform.t_in_years = waveform.t / 31556926;
waveform.a = SM(1,2);
waveform.e = SM(:,3);
waveform.p = SM(:,4);
waveform.iota_deg = SM(:,5);

% record the initial parameters
waveform.e0 = waveform.e(1);
waveform.p0 = waveform.p(1);
waveform.iota0_deg = waveform.iota_deg(1);
waveform.iota_hughes_deg0 = SM(1,12);
waveform.M = 1e6;
waveform.mu = 10;
[waveform.BigSteps ~] = size(SM);
waveform.SecPerMsun = 4.9255e-6;
waveform.SecPerM = waveform.SecPerMsun * waveform.M;
waveform.SecPermu = waveform.SecPerMsun * waveform.mu;
waveform.t_in_M = waveform.t / waveform.SecPerM; 
waveform.tol = 1e-4;

% fields that I'm not sure how I should set yet
%waveform.tol = 1e-2;
%waveform.t0 = 0;
%waveform.r0 = r0;
%waveform.theta0_deg = theta0_deg;
%waveform.phi0_deg = phi0_deg;
%waveform.sign_rdot0=sign_rdot0;
%waveform.sign_Tdot0=sign_Tdot0;
%[waveform.thetaE_deg waveform.phiE_deg] = EarthSkyPosition(S);
%if waveform.phiE_deg ~= 0
%    error('Your phi_E is not zero.  We are currently using Curt''s formula for the waveform, which requires phi_E = 0.');
%end

% for each time, we read inf a Teukolsky snapshot
for j = 1:waveform.BigSteps
    
    tic
    
    % status report
    display(['working on snapshot number ' num2str(j) ' / ' num2str(waveform.BigSteps)]);
    
    % calculate a KerrGeodesic structure for this orbit, even though it
    % isn't used to calculate H
    waveform.orbit{j} = KerrGeodesic(waveform.a, waveform.e(j), ...
        waveform.p(j), waveform.iota_deg(j), waveform.tol);
    %waveform.orbit{j} = GetFrequencies(waveform.a, waveform.e(j), ...
    %   waveform.p(j), waveform.iota_deg(j), waveform.tol);
    
    
    % build current data filename from SnapshotMatrix
    fname = sprintf('/data/matKITS/9.000000e-01_%e_%e_%e_1.000000e-02.mat',...
        waveform.e(j),waveform.p(j),waveform.iota_deg(j));
    load(fname);
    
    % compute the new entry to H that goes along with this data
    waveform.H{j} = HfromT(T,waveform.orbit{j}.OmegaR,...
        waveform.orbit{j}.OmegaTheta, waveform.orbit{j}.OmegaPhi, theta_rad);
    
    toc
    
end


% log the computational cost of this job
waveform.CPUsec = cputime - InitialCPUTime;

