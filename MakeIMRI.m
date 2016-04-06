%function MakeIMRI
%
%
%

% initialize structure
IMRIParameters

% initialize the waveform
waveform = InitializeWaveformSeries(S);

% spline it
output = SplineWaveformSeries(waveform);

% how do I observe it?
hin = 1;
r=1;
theta=45;
phi=0;
order='quadrupole';
[hplus hcross]=SeriesObserveWaveform(hin, r, theta, phi, order)