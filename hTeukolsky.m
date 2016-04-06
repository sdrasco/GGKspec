function h = hTeukolsky(filename,phi)
% h = hTeukolsky(filename)
% 
% Input is a file name pointing to data from Drasco & Hughes Teukolsky data
% with the format:
%
%  Col#  quantity
%
%   1       l
%   2       m
%   3       k
%   4       n
%   5       omega
%   6       Z^H.real()
%   7       ZH.imag()
%   8       Z^Inf.real()
%   9       ZInf.imag()
%   10       alpha 
%   11      S(pi/2-Teps) 
%   12      S(pi/3) 
%   13      S(pi/6) 
%   14      S(Teps) 
%   15      EInf 
%   16      EH 
%   17      LInf 
%   18      LH
%   19      Ain.real() 
%   20      Ain.imag()
%
% Output is the series object for
%
% h = (h_+) + i (h_x)
%
% Steve Drasco

% load data
T = load(filename);

% degenerate mkn matrix
S.mkn = T(:,2:4);

% strip out frequencies
omega = T(:,5);

% get the fundamental frequencies
[tf ind] = ismember([1 0 0],S.mkn,'rows');
if tf 
    OmegaPhi = omega(ind);
else
    error('failed to extract phi frequency');
end
[tf ind] = ismember([0 1 0],S.mkn,'rows');
if tf
    OmegaTheta = omega(ind);
else
    error('failed to extract theta frequency');
end
[tf ind] = ismember([1 0 1],S.mkn,'rows');
if tf
   OmegaR = omega(ind) - OmegaPhi;
else
    error('failed to extract radial frequency');
end
h.freqs = [OmegaPhi; OmegaTheta; OmegaR];

% Zlmkn
Zlmkn = T(:,6) + (1i) * T(:,7);

% Slmkn(pi/3);
Slmkn = T(:,12);

% m
m = S.mkn(:,1);

% hlmkn(pi/3)
S.coefs = (-2./(omega.^2)) .* exp((1i) * m * phi) .* Zlmkn .* Slmkn;

% complete degenerate structure
S.N = length(S.coefs);

% eliminate any terms with an absolute value less than 10 orders of
% magnitude less than the biggest absolute value
% abscoefs = abs(S.coefs);
% KillThreshold = 1e-10;
% kill = abscoefs < KillThreshold * max(abscoefs);
% S.coefs(kill) = [];
% S.mkn(kill,:) = [];
% S.N = length(S.coefs);

% make a list of the unique modes that will appear in s12
[h.mkn I J] = unique(S.mkn,'rows');
h.N = length(h.mkn(:,1));

% the vector J is such that degenerate_mkn = s12.mkn(J,:).  In this block we use
% this to collect degenerate modes
for L = 1:h.N
    ind = J == L;
    h.coefs(L) = sum(S.coefs(J == L));
end
h.coefs = h.coefs(:);

% final truncation
h = TruncateSeries(h,1e-10);
