function [tf s12kn] = SeriesProductElement(s1,s2,k,n)
%
% s12kn : k-n element of product of series s1 and s2
%
% Steve Drasco

% make sure these series are using the same frequencies, and copy
% fundamental frequencies into output structure
if ~isequal(s1.freqs, s2.freqs)
    error('SeriesProductElement: Disagreement in fundamental frequencies');
end

% find the modes in the product with the desired frequency
s12kn = 0;
for L = 1:s2.N
    kind = k == s2.mkn(L,2) + s1.mkn(:,2);
    nind = n == s2.mkn(L,3) + s1.mkn(:,3);
    ind = 2 == kind+nind;
    s12kn = s12kn + sum(s1.coefs(ind) * s2.coefs(L));
end

% return the result
if s12kn == 0
    tf = false;
else
    tf = true;
end
