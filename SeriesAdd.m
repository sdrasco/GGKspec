function s12 = SeriesAdd(s1,s2)
%
% The series s12 is the addition of series s1 and series s2.
% 
% Each series has fields:
%
%       N = number of harmonics in the series
%   freqs = 3 x 1 vector of fundamental frequencies
%     mkn = N x 3 matrix of indices (m, k, n)
%   coefs = N x 1 vector of coefficients 
%
% Such that
% 
%   s(x) = sum_{q=1}^N  s.coefs(q) exp(-i( (s.mkn(q,:) .* freqs) * x )
%
% Steve Drasco

% make these series have the same frequencies
if isequal(s1.freqs, s2.freqs)
    s12.freqs = s1.freqs;
else
    error('SeriesAdd: Disagreement in fundamental frequencies');
end

% add amplitudes of common modes
[s12.mkn ind1 ind2] = intersect(s1.mkn,s2.mkn,'rows');
s12.coefs = s1.coefs(ind1) + s2.coefs(ind2);

% strip common modes out of s1 and s2
s1.mkn(ind1,:) = [];
s1.coefs(ind1) =[];
s2.mkn(ind2,:) = [];
s2.coefs(ind2) = [];

% copy remaining modes into s12
s12.mkn = [s12.mkn; s1.mkn; s2.mkn];
s12.coefs = [s12.coefs; s1.coefs; s2.coefs];

% count the number of modes
s12.N = length(s12.coefs);

% sort summed series in order of increasing frequencies
s12 = SeriesSort(s12);

% truncate
s12 = TruncateSeries(s12,1e-10);

