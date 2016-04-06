function output = SplineWaveformSeries(waveform)
%
%
% Steve Drasco

% start clock to keep track of computational cost
InitialCPUTime = cputime;

% spline all the simple principle orbital element fields in the geodesics
e = zeros(waveform.BigSteps,1);
p = zeros(waveform.BigSteps,1);
iota_deg = zeros(waveform.BigSteps,1);
hughes_iota_deg = zeros(waveform.BigSteps,1);
E = zeros(waveform.BigSteps,1);
L = zeros(waveform.BigSteps,1);
Q = zeros(waveform.BigSteps,1);
rmin = zeros(waveform.BigSteps,1);
rmax = zeros(waveform.BigSteps,1);
theta_ = zeros(waveform.BigSteps,1);
OmegaPhi = zeros(waveform.BigSteps,1);
OmegaTheta = zeros(waveform.BigSteps,1);
OmegaR = zeros(waveform.BigSteps,1);
Gamma = zeros(waveform.BigSteps,1);
UpsilonPhi = zeros(waveform.BigSteps,1);
UpsilonTheta = zeros(waveform.BigSteps,1);
UpsilonR = zeros(waveform.BigSteps,1);
LP = zeros(waveform.BigSteps,1);
LT = zeros(waveform.BigSteps,1);
LR = zeros(waveform.BigSteps,1);
for i=1:waveform.BigSteps
    output.e(i) = waveform.orbit{i}.e;
    output.p(i) = waveform.orbit{i}.p;
    output.iota_deg(i) = waveform.orbit{i}.iota_deg;
    output.hughes_iota_deg(i) = waveform.orbit{i}.hughes_iota_deg;
    output.E(i) = waveform.orbit{i}.E;
    output.L(i) = waveform.orbit{i}.L;
    output.Q(i) = waveform.orbit{i}.Q;
    output.rmin(i) = waveform.orbit{i}.rmin;
    output.rmax(i) = waveform.orbit{i}.rmax;
    output.theta_(i) = waveform.orbit{i}.theta_;
    output.OmegaPhi(i) = waveform.orbit{i}.OmegaPhi;
    output.OmegaTheta(i) = waveform.orbit{i}.OmegaTheta;
    output.OmegaR(i) = waveform.orbit{i}.OmegaR;
    output.Gamma(i) = waveform.orbit{i}.Gamma;
    output.UpsilonPhi(i) = waveform.orbit{i}.UpsilonPhi;
    output.UpsilonTheta(i) = waveform.orbit{i}.UpsilonTheta;
    output.UpsilonR(i) = waveform.orbit{i}.UpsilonR;
    output.LP(i) = waveform.orbit{i}.LP;
    output.LT(i) = waveform.orbit{i}.LT;
    output.LR(i) = waveform.orbit{i}.LR;
end
output.e_pp = spline(waveform.t_in_M, output.e);
output.p_pp = spline(waveform.t_in_M, output.p);
output.iota_deg_pp = spline(waveform.t_in_M, output.iota_deg);
output.hughes_iota_deg_pp = spline(waveform.t_in_M, output.hughes_iota_deg);
output.E_pp = spline(waveform.t_in_M, output.E);
output.L_pp = spline(waveform.t_in_M, output.L);
output.Q_pp = spline(waveform.t_in_M, output.Q);
output.rmin_pp = spline(waveform.t_in_M, output.rmin);
output.rmax_pp = spline(waveform.t_in_M, output.rmax);
output.theta__pp = spline(waveform.t_in_M, output.theta_);
output.OmegaPhi_pp = spline(waveform.t_in_M, output.OmegaPhi);
output.OmegaTheta_pp = spline(waveform.t_in_M, output.OmegaTheta);
output.OmegaR_pp = spline(waveform.t_in_M, output.OmegaR);
output.Gamma_pp = spline(waveform.t_in_M, output.Gamma);
output.UpsilonPhi_pp = spline(waveform.t_in_M, output.UpsilonPhi);
output.UpsilonTheta_pp = spline(waveform.t_in_M, output.UpsilonTheta);
output.UpsilonR_pp = spline(waveform.t_in_M, output.UpsilonR);
output.LP_pp = spline(waveform.t_in_M, output.LP);
output.LT_pp = spline(waveform.t_in_M, output.LT);
output.LR_pp = spline(waveform.t_in_M, output.LR);


% now get a complete list of modes so that we can spline all of those too
mkn = waveform.H{1}.mkn;
for i=2:waveform.BigSteps
    mkn = union(mkn,waveform.H{i}.mkn,'rows');
end

% spline each mode, setting the coefficient value to zero for any time 
% that the mode was not present
[Nmodes ~] = size(mkn);
for mode_index=1:Nmodes 
    %display(['working on mode ' num2str(mode_index) ' out of ' num2str(Nmodes)]);
    coefs{mode_index} = zeros(waveform.BigSteps,1);
    for time_index=1:waveform.BigSteps
        % this block of logic calculations is the same as
        %    [~ ind] = ismember(mkn(mode_index,:),waveform.H{time_index}.mkn,'rows');
        % but is about twice as fast
        ind = mkn(mode_index,1) == waveform.H{time_index}.mkn(:,1);
        ind = ind + (mkn(mode_index,2) == waveform.H{time_index}.mkn(:,2));
        ind = ind + (mkn(mode_index,3) == waveform.H{time_index}.mkn(:,3));
        ind = ind == 3;
        if find(ind)
            coefs{mode_index}(time_index) = waveform.H{time_index}.coefs(ind);
        end
    end
    coefs_pp{mode_index} = spline(waveform.t_in_M, coefs{mode_index});
end


% organize the output
output.t = waveform.t_in_M;
output.mkn = mkn;
output.coefs = coefs;
output.coefs_pp = coefs_pp;
output.CPUsec = cputime - InitialCPUTime;