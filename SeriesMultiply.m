function s12 = SeriesMultiply(s1,s2,kmax,nmax)
%
% s12 : series structure for product of series s1 and s2
%
% each series structure has the following fields:
%
%       N = number of harmonics in the series
%   freqs = 3 x 1 vector of fundamental frequencies
%     mkn = N x 3 matrix of indices (m, k, n)
%   coefs = N x 1 vector of coefficients 
%
% Such that
% 
%   s1(t) = sum_{q=1}^N  s1.coefs(q) exp(-i( s1.mkn(q,:)*w )
%
% where
%
%   w = [w_phi; w_theta; w_r]
%
% Steve Drasco

% make sure these series are using the same frequencies, and copy
% fundamental frequencies into output structure
if isequal(s1.freqs, s2.freqs)
    s12.freqs = s1.freqs;
else
    error('SeriesMultiply: Disagreement in fundamental frequencies');
end

% get a degenerate list of harmonics in product of s1 and s2
degenerate_mkn = zeros(s2.N * s1.N,3);
degenerate_coefs = zeros(s2.N * s1.N,1);
for L = 1:s2.N
    RowRange = ((L-1)*s1.N + 1):L*s1.N;
    degenerate_mkn(RowRange,:) = [s2.mkn(L,1)+s1.mkn(:,1) s2.mkn(L,2)+s1.mkn(:,2) s2.mkn(L,3)+s1.mkn(:,3)];
    degenerate_coefs(RowRange) = s1.coefs * s2.coefs(L);
end

% eliminate any terms with an absolute value less than 10 orders of
% magnitude less than the biggest absolute value
abscoefs = abs(degenerate_coefs);
KillThreshold = 1e-10;
kill = abscoefs < KillThreshold * max(abscoefs);
degenerate_coefs(kill) = [];
degenerate_mkn(kill,:) = [];

% if we are restricting product series to a certain set of modes, eliminate
% uninteresting modes now
if nargin > 2
    KillIndex = find(abs(degenerate_mkn(:,2)) > kmax);
    degenerate_mkn(KillIndex,:) = [];
    degenerate_coefs(KillIndex) = [];
    KillIndex = find(abs(degenerate_mkn(:,3)) > nmax);
    degenerate_mkn(KillIndex,:) = [];
    degenerate_coefs(KillIndex) = [];
end

% make a list of the unique modes that will appear in s12
[s12.mkn I J] = unique(degenerate_mkn,'rows');
s12.N = length(s12.mkn(:,1));

% the vector J is such that degenerate_mkn = s12.mkn(J,:).  In this block we use
% this to collect degenerate modes
for L = 1:s12.N
    ind = J == L;
    s12.coefs(L) = sum(degenerate_coefs(J == L));
end
s12.coefs = s12.coefs(:);

% final truncation
s12 = TruncateSeries(s12,1e-6);
