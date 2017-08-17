%% Flux ratio vs koff at several affinities (koff and KD independent)
% For Fig. 2 of paper (flexible linker form of DB used)
% Does not incorporate koff2 (sliding koff)

% Re-dimensionalize without koff2
flux_3 = subs(subs(subs(subs(flux_2, DB, (D*alpha)/((alpha+3*D)*koff)),...
    alpha, ll*koff),gam, kon*Nt/koff), DF, D/(koff));
flux_3 = simplify(simplify(flux_3));

% No more numerical values need to be defined in this case
kon_n = 10^(-3);
% Plug in all numerical values and set kon = koff/KD
% x = ll; y = KD.
flux_final(x,y) = simplify(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, x),...
    kon, kon_n),koff,10^(-3)*y));

% Find the flux when there's no binding. (Should be a constant)
flux_no_binding = simplify(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    koff, koff),kon,0));

% Normalize the flux to the non-binding flux.
flux_ratio(x,y) = flux_final(x,y)/flux_no_binding;

close all

% Make an x-axis for the plot (do not start at zero; it gets mad).
x_axis = logspace(1,3);

close all
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 700, 600]);
% Plot flux ratio vs Dtilde/D for several KD values
for kd_n=[0.01, 0.1,1,10] % in uM
flux_ratio_final(x) = flux_ratio(x,kd_n); %#ok<SAGROW>
y_values = double(flux_ratio_final(x_axis));
semilogx(x_axis,y_values,'LineWidth',3);
xlabel('L_cl_p (nm^2)','fontsize',fontsize);
ylabel('Selectivity','fontsize',fontsize);
ax=gca;
ax.XTick = [10 100 1000];
ax.YTick = [10 20 30 40 50];
ax.XTickLabel = {[ fs ' 10'],[ fs ' 100'],[ fs ' 1000']};
ax.YTickLabel = {[ fs ' 10'],[ fs ' 20'],[ fs ' 30'],[ fs ' 40'],[ fs ' 50']};
hold all
end
% h = legend({'10 nM','100 nM', '1 uM', '10 uM'});
% title(h,'K_D')
% title('Normalized flux vs. L_cl_p')


%% Flux ratio vs koff at several affinities (koff and KD independent)
% For Fig. 2 of paper (flexible linker form of DB used)
% Does not incorporate koff2 (sliding koff)

% Re-dimensionalize without koff2
flux_3 = subs(subs(subs(subs(flux_2, DB, (D*alpha)/((alpha+3*D)*koff)),...
    alpha, ll*koff),gam, kon*Nt/koff), DF, D/(koff));
flux_3 = simplify(simplify(flux_3));

% No more numerical values need to be defined in this case
kon_n = 10^(-3);
% Plug in all numerical values and set kon = koff/KD
% x = ll; y = KD.
flux_final(x,y) = simplify(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, x),...
    kon, kon_n),koff,10^(-3)*y));

% Find the flux when there's no binding. (Should be a constant)
flux_no_binding = simplify(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    koff, koff),kon,0));

% Normalize the flux to the non-binding flux.
flux_ratio(x,y) = flux_final(x,y)/flux_no_binding;

close all

% Make an x-axis for the plot (do not start at zero; it gets mad).
x_axis = logspace(-2,2);

close all
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 700, 600]);
% Plot flux ratio vs Dtilde/D for several KD values
for lll_n=[10, 100,500,1000] % in uM
flux_ratio_final(y) = flux_ratio(lll_n,y); %#ok<SAGROW>
y_values = double(flux_ratio_final(x_axis));
semilogx(x_axis,y_values,'LineWidth',3);
xlabel('K_D (\muM)','fontsize',fontsize);
ylabel('Selectivity','fontsize',fontsize);
ax=gca;
ax.XTick = [0.01 0.1 1 10 100 1000];
ax.YTick = [10 20 30 40 50];
ax.XTickLabel = {[ fs ' 10^{-2}'],[ fs ' 10^{-1}'],[ fs ' 1'],[ fs ' 10'],[ fs ' 10^2'],[ fs ' 10^3']};
ax.YTickLabel = {[ fs ' 10'],[ fs ' 20'],[ fs ' 30'],[ fs ' 40'],[ fs ' 50']};
hold all
hold all
end
h = legend({'10','100', '500', '1000'},'fontsize',fontsize);
title(h,'L_cl_p (nm^2)','fontsize',fontsize)
% title('Normalized flux vs. K_D')


%% Flux ratio vs koff at several affinities (koff and KD independent)
% For Fig. 2 of paper (flexible linker form of DB used)
% Does not incorporate koff2 (sliding koff)
close all

% Make an x-axis for the plot (do not start at zero; it gets mad).
x_axis = linspace(0,500);

close all

% Plot flux ratio vs Dtilde/D for several KD values
for koff=10^(-3)*[0.1,1,10,100] % in uM
plot(x_axis,1.*koff.*x_axis.*1./(3.*1+ koff.*x_axis.*1),'LineWidth',3);
xlabel('L_c (nm)');
ylabel('D_B/D');
%axis([0 D_n 0 60])
%axis([0 0.001 0 10]); % use this to show that all curves start at 1
hold all
end
h = legend({'0.1', '1', '10','100'});
title(h,'K_D (\mu M)')
title('Diffusion coefficient ratio vs. peptide contour length')

%% Flux ratio vs koff at several affinities (koff and KD independent)
% For Fig. 2 of paper (flexible linker form of DB used)
% Does not incorporate koff2 (sliding koff)

% Re-dimensionalize without koff2
flux_3 = subs(subs(subs(subs(flux_2, DB, (D*alpha)/((alpha+3*D)*koff)),...
    alpha, ll*koff),gam, kon*Nt/koff), DF, D/(koff));
flux_3 = simplify(simplify(flux_3));

% No more numerical values need to be defined in this case
kon_n = 10^(-3);
% Plug in all numerical values and set kon = koff/KD
% x = ll; y = KD.
flux_final(x,y) = simplify(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, x),...
    kon, kon_n),koff,10^(-3)*y));

% Find the flux when there's no binding. (Should be a constant)
flux_no_binding = simplify(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    koff, koff),kon,0));

% Normalize the flux to the non-binding flux.
flux_ratio(x,y) = flux_final(x,y)/flux_no_binding;

close all

% Make an x-axis for the plot (do not start at zero; it gets mad).
x_axis = logspace(-2,2);

close all

% Plot flux ratio vs Dtilde/D for several KD values
for lll_n=[10, 100,500,1000] % in uM
flux_ratio_final(y) = flux_ratio(lll_n,1/y); %#ok<SAGROW>
y_values = double(flux_ratio_final(x_axis));
semilogx(x_axis,y_values,'LineWidth',3);
xlabel('K_A (\mu M^{-1})');
ylabel('Normalized flux');
%axis([0 D_n 0 60])
%axis([0 0.001 0 10]); % use this to show that all curves start at 1
hold all
end
h = legend({'10','100', '500', '1000'});
title(h,'L_c (nm)')
title('Normalized flux vs. association constant')
