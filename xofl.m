function x=xofl(orbit, lambda, t0, r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0)
%
% x=xofl(orbit, lambda, t0, r0, theta0_deg, phi0_deg, sign_rdot0, sign_Tdot0)
%
% This function is identical to KludgedXofL except that its input is only a
% single geodesic.  That is, the input called "inspiral" in the kludeged
% routine is replaced here with an "orbit" structure from KerrGeodesic.
%
% Returns Boyer-Lindquist (BL) coordinates 
%
%         x.t = BL t/M
%         x.r = BL r/M
%     x.theta = BL theta
%       x.phi = BL phi
%    x.CPUsec = number of seconds spent building this structure
%
% as a function of lambda, given the following inputs:
%
%  orbit: output structure from GEODESIC
% lambda: array of lambda values at which you want positions
%
% Output also has first, second, and third derivatives of the world 
% line.  Note that t is differentiated with respect to lambda, and 
% while spatial coordinates are differentiated with respect to t.
%
%         x.dt = dt/dlambda
%        x.ddt = d/dlambda(dt/dlambda)
%       x.dddt = d/dlambda(d/dlambda(dt/dlambda))
%         x.dr = dr/dt
%        x.ddr = d/dt(dr/dt)
%       x.dddr = d/dt(d/dt(dr/dt))
%     x.dtheta = dtheta/dt
%    x.ddtheta = d/dt(dtheta/dt)
%   x.dddtheta = d/dt(d/dt(dtheta/dt))
%       x.dphi = dphi/dt
%      x.ddphi = d/dt(dphi/dt)
%     x.dddphi = d/dt(d/dt(dphi/dt))
%
% See also KLUDGEDXOFL KERRGEODESIC 
%
% Steve Drasco
%

% start clock to keep track of computational cost
InitialCPUTime = cputime;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   translate initial coordinates to initial phases   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% set up options structure for fzero (tells how accuratly we'll get any of
% the phases should they be obtained via root finding).
tol=orbit.tol;
options = optimset('TolX',min([tol 1e-5]));

% get initial lambda_r (by root finding if not starting at turning point)
if abs(orbit.rmin - r0) < tol * orbit.rmin
    lambda_r = 0;
elseif abs(orbit.rmax - r0) < tol * orbit.rmax
    lambda_r = orbit.LR/2;
else
    nYr = (0:orbit.nmax)' * orbit.UpsilonR;
    rn = orbit.rn';
    rn(1) = rn(1)/2;
    fun_r =@(lambda_r) r0 - (2*rn*cos(nYr * lambda_r));
    if sign_rdot0 == 1 
        lambda_r = fzero(fun_r,orbit.LR*[0.5 1],options);
    else
        lambda_r = fzero(fun_r,orbit.LR*[0 0.5],options);
    end
end

% get initial lambda_theta (by root finding if not starting at turning
% point)
if abs(orbit.theta_ - (theta0_deg*pi/180)) < tol*orbit.theta_
    lambda_T = 0;
elseif abs((pi-orbit.theta_) - (theta0_deg*pi/180)) < tol*(pi-orbit.theta_)
    lambda_T = orbit.LT/2;
else
    kYT = (0:orbit.kmax)' * orbit.UpsilonTheta;
    thetak = orbit.thetak';
    thetak(1) = thetak(1)/2;
    fun_T =@(lambda_T) (theta0_deg*pi/180) - (2*thetak*cos(kYT * lambda_T));
    if sign_Tdot0 == 1 
        lambda_T = fzero(fun_T,orbit.LT*[0.5 1],options);
    else
        lambda_T = fzero(fun_T,orbit.LT*[0 0.5],options);
    end
end

% get initial lambda_phi (analytically)
nYr = (1:orbit.nmax)' * orbit.UpsilonR;
kYT = (1:orbit.kmax)' * orbit.UpsilonTheta;
phin = orbit.phin';
phik = orbit.phik';
lambda_phi = phi0_deg * pi/180;
lambda_phi = lambda_phi + phin*sin(nYr*lambda_r);
lambda_phi = lambda_phi + phik*sin(kYT*lambda_T);
lambda_phi = -lambda_phi / orbit.UpsilonPhi;

% get initial lambda_t (analytically)
tn = orbit.tn';
tk = orbit.tk';
lambda_t = t0;
lambda_t = lambda_t + tn*sin(nYr*lambda_r);
lambda_t = lambda_t + tk*sin(kYT*lambda_T);
lambda_t = -lambda_t / orbit.Gamma;

% put initial phases into a vector
lambda_x0 = [lambda_t lambda_r lambda_T lambda_phi];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Now we do the main calculation of the world line and its derivatives.   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% r and its contribution to t and phi %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% evolove positional element by requiring dw/dlambda = Upsilon
w0 = -orbit.UpsilonR * lambda_x0(2);
w = w0 + orbit.UpsilonR * lambda;

% compute exponential
Ex = exp(-(1i)*w);
clear w;

% loop over radial index
x.r = zeros(size(lambda));
dr = zeros(size(lambda));
ddr = zeros(size(lambda));
dddr = zeros(size(lambda));
phir = zeros(size(lambda));
dphir = zeros(size(lambda));
ddphir = zeros(size(lambda));
dddphir = zeros(size(lambda));
tr = zeros(size(lambda));
dtr = zeros(size(lambda));
ddtr = zeros(size(lambda));
dddtr = zeros(size(lambda));
for j=1:orbit.nmax;
    
    % compute cosine and sine terms
    cosjw = real(Ex.^j);
    sinjw = -imag(Ex.^j);
    
    % nth-contribution to t and its derivatives
    coef = orbit.tn(j);
    sinterm = coef .* sinjw;
    costerm = coef .* cosjw;
       tr = tr    + sinterm;
      dtr = dtr   + j*costerm;
     ddtr = ddtr  - (j^2 * sinterm);
    dddtr = dddtr - (j^3 * costerm);
    
    % nth-contribution to phi and its derivatives
    coef = orbit.phin(j);
    sinterm = coef .* sinjw;
    costerm = coef .* cosjw;
       phir = phir    + sinterm;
      dphir = dphir   + j*costerm;
     ddphir = ddphir  - (j^2 * sinterm);
    dddphir = dddphir - (j^3 * costerm);
    
    % nth-contribution to r and its derivatives
    % note that orbit.rn starts from n = 0
       coef = 2 * orbit.rn(j+1);
    sinterm = coef .* sinjw;
    costerm = coef .* cosjw;
      x.r = x.r  + costerm;
       dr = dr   - j*sinterm;  
      ddr = ddr  - (j^2)*costerm;  
     dddr = dddr + (j^3)*sinterm;  
end

% clear three long vectors to make room for frequency and its powers
clear coef sinterm costerm;

% scale by appropriate powers of Mino-frequencies
Y =  orbit.UpsilonR;
Y2 = Y.^2;
Y3 = Y.^3;
dtr = dtr .* Y;
ddtr = ddtr .* Y2;
dddtr = dddtr .* Y3;
dphir = dphir .* Y;
ddphir = ddphir .* Y2;
dddphir = dddphir .* Y3;
dr = dr .* Y;
ddr = ddr .* Y2;
dddr = dddr .* Y3;
clear Y Y2 Y3;

% add zeroth term to radius
x.r = x.r + orbit.rn(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theta and its contribution to t and phi %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% evolove positional element
w0 = -orbit.UpsilonTheta * lambda_x0(3);
w = w0 + orbit.UpsilonTheta*lambda;

% compute exponential
Ex = exp(-(1i)*w);
clear w;

% loop over radial index
x.theta = zeros(size(lambda));
dtheta = zeros(size(lambda));
ddtheta = zeros(size(lambda));
dddtheta = zeros(size(lambda));
phitheta = zeros(size(lambda));
dphitheta = zeros(size(lambda));
ddphitheta = zeros(size(lambda));
dddphitheta = zeros(size(lambda));
ttheta = zeros(size(lambda));
dttheta = zeros(size(lambda));
ddttheta = zeros(size(lambda));
dddttheta = zeros(size(lambda));
for j=1:orbit.nmax;
    
    % compute cosine and sine terms
    cosjw = real(Ex.^j);
    sinjw = -imag(Ex.^j);
    
    % nth-contribution to t and its derivatives
    coef = orbit.tk(j);
    sinterm = coef .* sinjw;
    costerm = coef .* cosjw;
       ttheta = ttheta    + sinterm;
      dttheta = dttheta   + j*costerm;
     ddttheta = ddttheta  - (j^2 * sinterm);
    dddttheta = dddttheta - (j^3 * costerm);
    
    % nth-contribution to phi and its derivatives
    coef = orbit.phik(j);
    sinterm = coef .* sinjw;
    costerm = coef .* cosjw;
       phitheta = phitheta    + sinterm;
      dphitheta = dphitheta   + j*costerm;
     ddphitheta = ddphitheta  - (j^2 * sinterm);
    dddphitheta = dddphitheta - (j^3 * costerm);
    
    % nth-contribution to r and its derivatives
    % note that orbit.rn starts from n = 0
       coef = 2 * orbit.thetak(j+1);
    sinterm = coef .* sinjw;
    costerm = coef .* cosjw;
      x.theta = x.theta  + costerm;
       dtheta = dtheta   - j*sinterm;  
      ddtheta = ddtheta  - (j^2)*costerm;  
     dddtheta = dddtheta + (j^3)*sinterm;  
end

% clear things that we're done with
clear Ex cosjw sinjw coef sinterm costerm;

% scale by appropriate powers of Mino-frequencies
Y =  orbit.UpsilonTheta;
Y2 = Y.^2;
Y3 = Y.^3;
dttheta = dttheta .* Y;
ddttheta = ddttheta .* Y2;
dddttheta = dddttheta .* Y3;
dphitheta = dphitheta .* Y;
ddphitheta = ddphitheta .* Y2;
dddphitheta = dddphitheta .* Y3;
dtheta = dtheta .* Y;
ddtheta = ddtheta .* Y2;
dddtheta = dddtheta .* Y3;
clear Y Y2 Y3;

% add zeroth term to theta 
% could also just do: x.theta = x.theta + pi/2
x.theta = x.theta + orbit.thetak(1);

%%%%%%%%%
%  phi  %
%%%%%%%%%

% evolve positional element
w0 = -orbit.UpsilonPhi * lambda_x0(4);
w = w0 + orbit.UpsilonPhi * lambda;

% add up remaining parts
x.phi = w + phitheta + phir;
clear phitheta phir;
dphi = orbit.UpsilonPhi + dphitheta + dphir;
clear dphitheta dphir;
ddphi = ddphitheta + ddphir;
clear ddphitheta ddphir;
dddphi = dddphitheta + dddphir;
clear dddphitheta dddphir;

%%%%%
% t %
%%%%%

% evolve positional element
w0 = -orbit.Gamma * lambda_x0(1);
w = w0 + orbit.Gamma * lambda;

% add up remaining parts
x.t = w + ttheta + tr;
clear ttheta tr;
x.dt = orbit.Gamma + dttheta + dtr;
clear dttheta dtr;
x.ddt = ddttheta + ddtr;
clear ddttheta ddtr;
x.dddt = dddttheta + dddtr;
clear dddttheta dddtr;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% final block to convert from lambda-derivatives to t-derivatives %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% We will need first five inverse powers of dt repeatedly. 
dtm1 = x.dt.^(-1);
dtm2 = x.dt.^(-2);
dtm3 = x.dt.^(-3);
dtm4 = x.dt.^(-4);
dtm5 = x.dt.^(-5);

% main calculation, freeing space as we go
x.dr = dr .* dtm1;
x.ddr = (ddr .* dtm2) - (dr .* x.ddt .* dtm3);
x.dddr = (dddr .* dtm3) - (3 * ddr .* x.ddt .* dtm4) ...
    - (dr .* x.dddt .* dtm4) + (3 * dr .* x.ddt.^2 .* dtm5);
clear dr ddr dddr;
x.dtheta = dtheta .* dtm1;
x.ddtheta = (ddtheta .* dtm2) - (dtheta .* x.ddt .* dtm3);
x.dddtheta = (dddtheta .* dtm3) - (3 * ddtheta .* x.ddt .* dtm4) ...
    - (dtheta .* x.dddt .* dtm4) + (3 * dtheta .* x.ddt.^2 .* dtm5);
clear dtheta ddtheta dddtheta;
x.dphi = dphi .* dtm1;
x.ddphi = (ddphi .* dtm2) - (dphi .* x.ddt .* dtm3);
x.dddphi = (dddphi .* dtm3) - (3 * ddphi .* x.ddt .* dtm4) ...
    - (dphi .* x.dddt .* dtm4) + (3 * dphi .* x.ddt.^2 .* dtm5);
clear dphi ddphi dddphi;

% Clear things just in case Matlab is lazy about it, since they're big
clear dtm1 dtm2 dtm3 dtm4

% store a copy of lambda in the output structure
x.lambda = lambda(:);

% log the computational cost of this job
x.CPUsec = cputime - InitialCPUTime;
