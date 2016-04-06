function [hplus, hcross, orbit] = SnapshotSpectrum(a, e, p, iota_deg, theta_obs_deg)

orbit = KerrGeodesic(a, e, p, iota_deg, 1e-10);
snapshot = SnapshotSeries(orbit);
H = ObserveSnapshot(snapshot,theta_obs_deg);
hplus = SeriesRe(H);
hcross = SeriesIm(H);
hcross.coefs = -hcross.coefs;

% plot the spectrum with lines
% figure
% hold on;
% S = hplus;
% zfloor = min(abs(S.coefs))/10;
% for i=1:length(S.coefs);
%     x = S.mkn(i,2);
%     y = S.mkn(i,3);
%     z = abs(S.coefs(i));
%     plot3([x x],[y y],[zfloor z],'b-');
% end
% axis([min(S.mkn(:,2)) max(S.mkn(:,2)) min(S.mkn(:,3)) max(S.mkn(:,3)) zfloor max(abs(S.coefs))]);
% set(gca,'zscale','log','fontsize',16);
% xlabel('k');
% ylabel('n'); 
% zlabel('| h_+^{mkn} |')
% box on;

% plot the spectrum with points
figure
S = hplus;
plot3(S.mkn(:,2),S.mkn(:,3),abs(S.coefs),'b.','markersize',10)
set(gca,'zscale','log','fontsize',16);
xlabel('k');
ylabel('n'); 
title('| h_+^{mkn} |')
box on;