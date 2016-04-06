function ImS = SeriesIm(S)
%
% The series ImS is the series representation of the imaginary part of 
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

% get negative conjugate S
ncS = SeriesConj(S);
ncS.coefs = -1*ncS.coefs;

% imaginary part is the (-i/2) * (S - conjS)
ImS = SeriesAdd(ncS,S);
ImS.coefs = -(1i)*ImS.coefs/2;