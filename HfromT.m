function H = HfromT(T, OmegaR, OmegaTheta, OmegaPhi, theta_rad)
%
% theta_rad must be one of: pi/2, pi/3, pi/6, 0
%

% get spin-weighted spherical harmonics
if theta_rad == pi/2
    S = T.S(:,1);
elseif theta_rad == pi/3
    S = T.S(:,2);
elseif theta_rad == pi/6
    S = T.S(:,3);
elseif theta_rad == 0
    S = T.S(:,4);
else
    error('theta_rad must be one of: pi/2, pi/3, pi/6, 0');
end

% set unique list of modes, and all other fields of output other than coefs
[H.mkn I J] = unique(T.lmkn(:,2:4),'rows');
[H.N ~] = size(H.mkn);
H.freqs = [OmegaPhi; OmegaTheta; OmegaR];

% compute degenerate list of coefs
omega = T.lmkn(:,2:4) * H.freqs;
Z = T.ZH;
degenerate_coefs = -2 * Z .* S ./ (omega.^2);

% the vector J is such that T.lmkn = H.mkn(J,:).  In this block we use
% this to collect degenerate modes
for L = 1:H.N
    ind = J == L;
    H.coefs(L) = sum(degenerate_coefs(J == L));
end
H.coefs = H.coefs(:);

% TEST code, leaving mode list degenerate
% H.mkn = T.lmkn(:,2:4);
% [H.N ~] = size(H.mkn);
% H.freqs = [OmegaPhi; OmegaTheta; OmegaR]; 
% omega = T.lmkn(:,2:4) * H.freqs;
% Z = T.ZH;
% H.coefs = -2 * Z .* S ./ (omega.^2);
% H.coefs = H.coefs(:);