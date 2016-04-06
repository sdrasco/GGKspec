function sofx = EvolvingSeriesEval(s,x)
%
% s is a structure representing an "evolving series" that
% incorporates coefficients and frequencies that change over time. 
% Structure s includes the following fields:
%
%          N = number of harmonics in the series
%        mkn = N x 3 matrix of indices (m, k, n)
%   coefs_pp = 1 x N array of spline structures of coefs 
%   freq1_pp = Spline structure corresponding to fundamental freq #1
%   freq2_pp = Spline structure corresponding to fundamental freq #2
%   freq3_pp = Spline structure corresponding to fundamental freq #3
%          t = Vector of times at which the splines are valid 
%              (eg. s.t = output.t) -- Necessary to ensure splines are only used
%              for interpolation.
%   
%
% Vector x corresponds to the values at which the series will be evaluated.
% When using tt, make sure tt is in the same units as output.t


if min(x) < min(s.t)
    error('EvolvingSeriesEval: Outside minimum time boundary. Time vector used to evaluate series must be within spline interpolation range.')
end

if max(x) > max(s.t)
    error('EvolvingSeriesEval: Outside maximum time boundary.  Time vector used to evaluate series must be within spline interpolation range.')
end

lengthx = length(x);

if s.N > 1000 && lengthx > 13000
    error('Number of time steps is too large for this number of coefficients. Memory usage issues will occur.')
end
 

% Spline coefs of each mode
for i = 1:s.N
    coefsSplined{i} = ppval(s.coefs_pp{i},x);
end

%%% Calculate exponent of series %%%

% Integrate frequency splines
intfreq1 = fnint(s.freq1_pp);
intfreq2 = fnint(s.freq2_pp);
intfreq3 = fnint(s.freq3_pp);

% Evaluate integrated freq splines at values of x
intfreq1Splined = ppval(intfreq1,x);
intfreq2Splined = ppval(intfreq2,x);
intfreq3Splined = ppval(intfreq3,x);

% Create frequency matrix [3 x length(x)]
intfreqs(1,:) = intfreq1Splined;
intfreqs(2,:) = intfreq2Splined;
intfreqs(3,:) = intfreq3Splined;

% Matrix multiply mkn matrix and frequency matrix to produce matrix of
% size [N x length(x)] representing each phi value for each mode (rows)
% at every time step (columns).
phi = s.mkn * intfreqs;

value2 = exp((1i)*phi);
   

% Calculate final vector by summing over all modes for each time value.
sofx = linspace(0,0,lengthx);
for i = 1:lengthx
    for j = 1:s.N
        totterm = coefsSplined{j}(i) * value2(j,i);
        sofx(i) = sofx(i) + totterm;
    end
end

end