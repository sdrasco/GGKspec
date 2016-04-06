function sout = SeriesSort(sin)
%
% Sorts the series modes and frequencies in order of increasing frequency
%
% Steve Drasco

% master matrix for sorting.  Note that we have to split off the real and
% imaginary parts individually because otherwise sortrows will sort in
% order of increasing absolute value
w = sin.mkn * sin.freqs;
M = [sin.mkn real(sin.coefs) imag(sin.coefs) w];

% sort the master matrix rows in order of incresing fifth-column-elements
M = sortrows(M,6);

% copy master matrix back into series structure
sout.mkn = M(:,1:3);
sout.coefs = M(:,4) + (1i)*M(:,5);
sout.freqs = sin.freqs;
sout.N = sin.N;

