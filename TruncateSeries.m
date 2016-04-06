function Sout = TruncateSeries(S,tol)
%
% simple routine to delete elements from a series for which the ratio of
% their magnitude to the term with the largest magnitude is less than a
% specified tollerance.
%

ind = abs(S.coefs) < tol * max(abs(S.coefs));
S.coefs(ind) = [];
S.mkn(ind,:) = [];
S.N = length(S.coefs);
Sout = S;

