function waveform = RemakeKludgeKITS()
% 
% 
%  Steve Drasco
% 

% start clock to keep track of computational cost
InitialCPUTime = cputime;

% load the Teukolsky inspiral structure
load KITS_TeukComplete_1e-3.mat

% undersample it 
NewBigSteps = 20;
index = round(linspace(1,waveform.BigSteps,NewBigSteps));
waveform.BigSteps = NewBigSteps;
waveform.t = waveform.t(index);
waveform.t_in_M = waveform.t_in_M(index);
waveform.t_in_years = waveform.t_in_years(index);
waveform.e = waveform.e(index);
waveform.p = waveform.p(index);
waveform.iota_deg = waveform.iota_deg(index);
OldOrbits = waveform.orbit;
waveform.thetaE_deg = 60;
rmfield(waveform,'orbit');

% for each time, we generate a snapshot
for j = 1:waveform.BigSteps
    RefIndex = index(j);
    tic
    display(['working on geodesic number ' num2str(j) ' / ' num2str(waveform.BigSteps)]);
    %waveform.orbit{j} = OldOrbits{RefIndex};
        waveform.orbit{j} = KerrGeodesic(waveform.a, waveform.e(j), ...
        waveform.p(j), waveform.iota_deg(j), waveform.tol);
    snapshot = SnapshotSeries(waveform.orbit{j});
    waveform.H{j} = ObserveSnapshot(snapshot,waveform.thetaE_deg);
    waveform.H{j}.coefs = waveform.H{j}.coefs; % * (waveform.mu * waveform.SecPerMsun) / S.D;
    toc
end

% log the computational cost of this job
waveform.CPUsec = cputime - InitialCPUTime;

