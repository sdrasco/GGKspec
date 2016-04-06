function ds = SeriesDifferentiate(s)
%
% ds is the derivative of the series s with respect to the time-argument of
% the series s.
%
% Steve Drasco

% copy s into ds
ds = s;

% get frequencies
w = ds.mkn * ds.freqs;

% rescale the coefficients in the derivative series
for L = 1:ds.N
    ds.coefs(L) = -(1i)*w(L)*s.coefs(L);
end

% truncate
ds = TruncateSeries(ds,1e-10);

