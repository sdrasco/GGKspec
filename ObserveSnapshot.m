function H = ObserveSnapshot(snapshot,theta_deg)
%
% There is a better routine in Serieshphx
%
% Steve Drasco
%

% convert angles to radians
theta = theta_deg * pi/180;

% term 1
term1 = snapshot.term1;
term1.coefs = (sin(theta))^2 * term1.coefs;

% term 2 
term2 = snapshot.term2;
term2.coefs = -1 * (1 + cos(theta)) * sin(theta) * term2.coefs;

% term 3 
term3 = snapshot.term3;
term3.coefs = (1 - cos(theta)) * sin(theta) * term3.coefs;

% term 4
term4 = snapshot.term4;
term4.coefs = 0.5 * (1 + (cos(theta))^2) * term4.coefs;

% term 5
term5 = snapshot.term5;
term5.coefs = cos(theta) * term5.coefs;

% combine terms
H = SeriesAdd(term1,term2);
H = SeriesAdd(H,term3);
H = SeriesAdd(H,term4);
H = SeriesAdd(H,term5);