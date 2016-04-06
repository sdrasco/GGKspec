function [h_mkn w t]=StitchMode(S,m,k,n,N)
%
% [h t]=StitchMode(S,m,k,n,N)
%
% Input structure is of the format returned by SplineWaveformSeries.
% The output is a time vector t made of N samples over the timespan S.t,
% and also the "single mode" of h_+ - i h_x
%
%   h = (h_+ - ih_x)_{mkn} = h_{mkn} * exp(-i w_{mkn})
% 
% where 
%
%   dw_{mkn}/dt = omega_{mkn}
%
% sampled N-times over the timespan of the waveform.
%
% Steve Drasco
%

% make time vector
t = linspace(S.t(1), S.t(end),N);

% For now, we take three initial phases to arbitrary values.  
% EDIT HERE ONCE I REMEMBER THE RELATIONSHIP BETWEEN lambda-phases and t-phases 
%
% Note: I think the relation is t_x = Gamma (lambda_x - lambda_t)
%
OmegaPhi0 = ppval(S.OmegaPhi_pp,t(1));
OmegaTheta0 = ppval(S.OmegaTheta_pp,t(1));
OmegaR0 = ppval(S.OmegaR_pp,t(1));
const1 = 0.387;
const2 = 0.813;
const3 = 0.50;
%tr0 = rand(1)*2*pi/OmegaPhi0;
%ttheta0 = rand(1)*2*pi/OmegaTheta0;
%tphi0 = rand(1)*2*pi/OmegaR0;
%Curt changed the random nums to fixed arbitrary constants, so theyd be
%same for every run
tr0 = const1*2*pi/OmegaPhi0;
ttheta0 = const2*2*pi/OmegaTheta0;
tphi0 = const3*2*pi/OmegaR0;

% evolve the part of the mode that depends only on the principle elements
[tf ind] = ismember([m k n],S.mkn,'rows');
if ~tf
    error('The mode you requested is not present in this waveform.');
end
h_mkn = ppval(S.coefs_pp{ind},t);

% get initial value of w_r, w_theta, and w_phi
wphi0 = -OmegaPhi0 * tphi0;
wtheta0 = -OmegaTheta0 * ttheta0;
wr0 = -OmegaR0 * tr0;

% get a spline representation of the integral the frequencies
IOmegaPhi_pp = fnint(S.OmegaPhi_pp);
IOmegaTheta_pp = fnint(S.OmegaTheta_pp);
IOmegaR_pp = fnint(S.OmegaR_pp);

% evolove positional elements by requiring dw/dlambda = Omega
wphi = wphi0 + (ppval(IOmegaPhi_pp, t) - ppval(IOmegaPhi_pp,0));
wtheta = wtheta0 + (ppval(IOmegaTheta_pp, t) - ppval(IOmegaTheta_pp,0));
wr = wr0 + (ppval(IOmegaR_pp, t) - ppval(IOmegaR_pp,0));
w = m*wphi + k*wtheta + n*wr;
clear wphi wtheta wr;

% return the result
%h = h_mkn .* exp(-(1i)*w); %Steve returned [h t], but Curt modified
%this to return amplitude and phase instead.



