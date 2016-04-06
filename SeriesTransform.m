function S = SeriesTransform(fseries,dtseries,eikwTDeltat,einwrDeltat,Gamma,tol,OnlyEvenk)
%
%
%
% Steve Drasco
%

if nargin == 6
    OnlyEvenk = false;
end

% get the product of dt and the series to be transformed
dtfseries = SeriesMultiply(dtseries,fseries);
dtfseries = TruncateSeries(dtfseries,tol);

% set index ranges leaving out zeros for special treatment
nrange = einwrDeltat{end};
krange = eikwTDeltat{end};
nrange(nrange==0) = [];
krange(krange==0) = [];

% if flagged to do so, we will only do the even k-modes
if OnlyEvenk
    krange = krange(mod(krange,2)==0);
end

% keep the zero-zero term if there is one
[tf,ind] = ismember([0 0 0], dtfseries.mkn,'rows');
if tf
    coefs = dtfseries.coefs(ind)/Gamma;
    mkn = [0 0 0];
else
    coefs = [];
    mkn = [];
end

% get the terms where both k and n are nonzero
for n=nrange
    nind = einwrDeltat{end}==n;
    nterm = SeriesMultiply(dtfseries,einwrDeltat{nind});
    nterm = TruncateSeries(nterm,tol);
    for k=krange
        kind = eikwTDeltat{end}==k;
        [tf Fkn] = SeriesProductElement(nterm,eikwTDeltat{kind},k,n);
        if tf
            coefs = [coefs; Fkn/Gamma];
            mkn = [mkn; 0 k n];
        end
    end
end

% get the terms where k is zero and n isn't
for n=nrange
    nind = einwrDeltat{end}==n;
    [tf Fkn] = SeriesProductElement(einwrDeltat{nind},dtfseries,0,n);
    if tf
        coefs = [coefs; Fkn/Gamma];
        mkn = [mkn; 0 0 n];
    end
end

% get the terms where n is zero and k isn't
for k=krange
    kind = eikwTDeltat{end}==k;
    [tf Fkn] = SeriesProductElement(eikwTDeltat{kind},dtfseries,k,0);
    if tf
        coefs = [coefs; Fkn/Gamma];
        mkn = [mkn; 0 k 0];
    end
end

% collect results
S.freqs = fseries.freqs/Gamma;
S.mkn = mkn;
S.coefs = coefs;
S.N = length(S.coefs);

% truncate
S = TruncateSeries(S,1e-10);
