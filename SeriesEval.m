function sofx = SeriesEval(s,x)
%
% the series s is a structure with the following fields:
%
%       N = number of harmonics in the series
%   freqs = 3 x 1 vector of fundamental frequencies
%     mkn = N x 3 matrix of indices (m, k, n)
%   coefs = N x 1 vector of coefficients 
%
% Such that
% 
%   sofx = sum_{q=1}^N  s.coefs(q) exp(-i( (s.mkn(q,:) .* freqs) * x )
%
%
% Steve Drasco

% frequencies
w = s.mkn * s.freqs;

% main loop over harmonics
sofx = zeros(size(x));
for L = 1:s.N
    
    % the new term 
    term = s.coefs(L) * exp(-(1i)*x*w(L));
        
    % increment sum
    sofx = sofx + term;
    
end