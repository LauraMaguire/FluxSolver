function results = fluxPrep()

disp('Loading files...');
r = LoadResults(0);
n = length(r.d);

disp('Calculating flux...');
for i=1:n
    [r.bindFlux(i), r.nonbindFlux(i)]= subNum(r.params{i},0);
    disp([num2str(i) ' of ' num2str(n) ' complete.']);
end

results = r;

end