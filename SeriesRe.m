function ReS = SeriesRe(S)
%
% The series ReS is the series representation of the real part of 
% the series S.
% 
% Each series has fields:
%
%       N = number of harmonics in the series
%   freqs = 3 x 1 vector of fundamental frequencies
%     mkn = N x 3 matrix of indices (m, k, n)
%   coefs = N x 1 vector of coefficients 
%
%
% Steve Drasco

% get conjugate S
cS = SeriesConj(S);

% real part is the (1/2) * (S + conjS)
ReS = SeriesAdd(cS,S);
ReS.coefs = ReS.coefs/2;

