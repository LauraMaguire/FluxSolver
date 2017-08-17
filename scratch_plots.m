figure('DefaultAxesFontSize',18)
hold all

hopValues = unique(r4.khop);
hop = cell(1,length(hopValues));
for i=1:length(hopValues)
    hop{i} = find(r4.khop == hopValues(i));
    leg{i} = num2str(hopValues(i));
    semilogx(r4.kd(hop{i}),r4.bindFlux(hop{i})./r4.nonbindFlux(hop{i}))
end
l = legend(leg,'Location','northeast');
v = get(l,'title');
set(v,'String','kHop')
xlabel('K_D (uM)');
ylabel('Selectivity');