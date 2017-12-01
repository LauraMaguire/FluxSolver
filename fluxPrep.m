function results = fluxPrep(r)

if nargin == 0
    disp('Loading files...');
    r = LoadResults(0);
end

n = length(r.deltat);

disp('Calculating flux...');
for i=1:n
    [r.bindFlux(i), r.nonbindFlux(i)]= subNum(r.params{i},0);
    disp([num2str(i) ' of ' num2str(n) ' complete.']);
    ptemp = r.params{i};
    ptemp.DB = r.dBound(i) + r.dErr(i);
    [r.bFluxErr(i), r.nFluxErr(i)] = subNum(ptemp,0);
end

results = r;

end