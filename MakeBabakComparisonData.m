function MakeBabakComparisonData()
%
%
%

% get the range in t used by stas et al
M=load('../BabakEtAlCheck/StasData1.dat');
t = M(:,1);
clear M;

% retrograde data set
% a = 0.9;
% e = 0.3;
% p = 12;
% iota_deg = 140;
% theta_deg = 60;
% phi = 0;
% tol = 1e-10;
% orbit = KerrGeodesic(a, e, p, iota_deg, tol);
% snapshot = SnapshotSeries(orbit);
% H = ObserveSnapshot(snapshot,theta_deg);
% k1 = [t real(SeriesEval(H,t))];
% save ../BabakEtAlCheck/QuadRetroDC.mat k1
% save Hretro.mat H orbit snapshot
% 
% display('finished retrograde');


% prograde data set
a = 0.9;
e = 0.7;
p = 6;
iota_deg = 60;
theta_deg = 90;
phi = 0;
tol = 1e-6;
orbit = KerrGeodesic(a, e, p, iota_deg, tol);
snapshot = SnapshotSeries(orbit);
H = ObserveSnapshot(snapshot,theta_deg);
k2 = [t real(SeriesEval(H,t))];
save ../BabakEtAlCheck/QuadProDC.mat k2
save Hpro.mat H orbit snapshot

% now the octupole
% 
% % retrograde data set
% a = 0.9;
% e = 0.3;
% p = 12;
% iota_deg = 140;
% theta_deg = 60;
% phi = 0;
% tol = 1e-10;
% orbit = KerrGeodesic(a, e, p, iota_deg, tol);
% octupole_flag = 1;
% h=SeriesCartesianMetric(orbit,octupole_flag);
% [hplus hcross]=SeriesObserveWaveform(h, 1, theta_deg, phi, 'octupole');
% k1 = [t real(SeriesEval(hplus,t))];
% save ../BabakEtAlCheck/RetroDC.mat k1
% 
% 
% % prograde data set
% a = 0.9;
% e = 0.7;
% p = 6;
% iota_deg = 60;
% theta_deg = 90;
% phi = 0;
% tol = 1e-4;
% display('starting KerrGeodesic')
% orbit = KerrGeodesic(a, e, p, iota_deg, tol);
% display('starting metric calculation');
% octupole_flag = 1;
% h=SeriesCartesianMetric(orbit,octupole_flag);
% display('starting sky projection');
% [hplus hcross]=SeriesObserveWaveform(h, 1, theta_deg, phi, 'octupole');
% hplus = TruncateSeries(hplus,1e-6);
% display('evaluating series');
% k2 = [t real(SeriesEval(hplus,t))];
% save ../BabakEtAlCheck/ProDC.mat k2
