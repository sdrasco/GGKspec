function conjS = SeriesConj(S)
%
% The series conjs is the complex conjugate of series s.
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


% copy frequencies and N
conjS.N = S.N;
conjS.freqs = S.freqs;

% switch sign of indices and conjugate coefs
conjS.mkn = -S.mkn;
conjS.coefs = conj(S.coefs);




