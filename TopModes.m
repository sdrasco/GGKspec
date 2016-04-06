function TopIndices = TopModes(S,num)
%
% S = output structure from SplineWaveformSeries
%

% get an array of the largest amplitude reached during the inspiral by each mode
N = length(S.coefs);
v = zeros(N,1);
for i=1:N
    v(i) = max(abs(S.coefs{i}));
end

% cycle through the largest amplitude vector getting the largest one and
% then zeroing it out before asking for the largest one again
TopTenIndices = zeros(num,1);
for i=1:num
    maxind = (v == max(v));
    new = find(maxind);
    new = new(1);
    TopIndices(i) =  new;
    v(maxind) = 0;
end
