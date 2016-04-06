function Hmkn = GetMode(waveform,mkn)
%
% Inputs are:
%   waveform: the structure returned by InitializeWaveformSeries
%        mkn: [m k n] mode that you want returned
%
% output is h_{mkn} for each snapshot contained in the waveform structure,
% given in the same order as waveform.orbit{}.
%
% Steve Drasco

% initialize
Hmkn=[];

% main loop to strip mode from each geodesic
for i=1:waveform.BigSteps
    [tf ind] = ismember(mkn,waveform.H{i}.mkn,'rows');
    omega = mkn*waveform.H{i}.freqs;
    if tf
        Hmkn = [Hmkn; waveform.H{i}.coefs(ind)];
    else
        Hmkn = [Hmkn; 0];
    end
end