%% Flux ratio vs kon for fixed DB, at several KDs
% For Fig. 1 of paper (no functional form of Dtilde)

% Re-dimensionalize without functional form of DB
flux_3 = subs(subs(subs(flux_2, DB, (DB/(koff))),...
    gam, kon*Nt/koff), DF, D/(koff));
flux_3 = simplify(simplify(flux_3));

%Define remaining numerical values.
kon_n = 10^(-3); % diffusion-limited on-rate

% Plug in all numerical values. DB = 0.1*D.
% x = Nt*Kon; y = KD.
flux_final(x,y) = simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
    flux_3,koff,kon*y ),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    DB, 0.1*D_n),Nt*kon, x),1,1));

flux_final = simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    DB, 5.88),koff,0.1),kon,kon_n));

% Find the flux when there's no binding. (Should be a constant)
flux_no_binding = simplify(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    koff, koff),kon,0));

% Normalize the flux to the non-binding flux.
flux_ratio(x,y) = flux_final(x,y)/flux_no_binding;

% Make an x-axis for the plot (do not start at zero; it gets mad).
x_axis = logspace(-4,1);

close all

% Plot flux ratio vs Dtilde/D for several KD values
for kd_n=[0.01, 0.1,1,10] % in uM
flux_ratio_final(x) = flux_ratio(x,kd_n); %#ok<SAGROW>
y_values = double(flux_ratio_final(x_axis));
semilogx(x_axis,y_values,'LineWidth',3);
xlabel('k_{on} (\mu s^{-1} \mu M^{-1})');
ylabel('Normalized flux');
%axis([0 D_n 0 60])
%axis([0 0.001 0 10]); % use this to show that all curves start at 1
hold all
end
h = legend({'10 nM','100 nM', '1 uM', '10 uM'});
title(h,'K_D')
title('Normalized flux vs. k_{on}')

%% Flux ratio vs Nt for fixed DB, at several KDs
% For Fig. 1 of paper (no functional form of Dtilde)
fontsize = 18;
fs = 18;
% Re-dimensionalize without functional form of DB
flux_3 = subs(subs(subs(flux_2, DB, (DB/(koff))),...
    gam, kon*Nt/koff), DF, D/(koff));
flux_3 = simplify(simplify(flux_3));

%Define remaining numerical values.
kon_n = 10^(-3); % diffusion-limited on-rate

% Plug in all numerical values. DB = 0.1*D.
% x = Nt; y = KD.
flux_final(x,y) = D_n*simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, x),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    DB, 0.1*D_n),koff,kon*y),kon,kon_n));

% flux_final = simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
%     flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
%     DB, 5.88),koff,0.1),kon,kon_n));

% Find the flux when there's no binding. (Should be a constant)
flux_no_binding = D_n*simplify(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    koff, koff),kon,0));

% Normalize the flux to the non-binding flux.
flux_ratio(x,y) = flux_final(x,y)/flux_no_binding;

% Make an x-axis for the plot (do not start at zero; it gets mad).
x_axis = logspace(1,4);

close all

% Plot flux ratio vs Dtilde/D for several KD values
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 700, 600]);
%figure('DefaultAxesFontSize',18)
for kd_n=[0.01, 0.1,1,10] % in uM
flux_ratio_final(x) = flux_ratio(x,kd_n); %#ok<SAGROW>
y_values = double(flux_ratio_final(x_axis));
semilogx(x_axis,y_values,'LineWidth',3);
xlabel('N_t (\muM)','fontsize',fontsize);
ylabel('Selectivity','fontsize',fontsize);
% axis([0 0.01 0 5])
% ax=gca;
% ax.XTick = [10 100 10^3 10^4];
% ax.YTick = [0 50 100 150];
% ax.XTickLabel = {[ fs ' 10'],[ fs ' 100'],[ fs ' 10^3'],[ fs ' 10^4']};
% ax.YTickLabel = {[ fs ' 0'],[ fs ' 50'],[ fs ' 100'],[ fs ' 150']};
hold all
end
%h = legend({'10 nM','100 nM', '1 uM', '10 uM'});
%title(h,'K_D')
%title('Normalized flux vs. N_t')

%% Flux ratio vs 1/D for fixed DB, at several KDs
% For Fig. 1 of paper (no functional form of Dtilde)

% Re-dimensionalize without functional form of DB
flux_3 = subs(subs(subs(flux_2, DB, (DB/(koff))),...
    gam, kon*Nt/koff), DF, D/(koff));
flux_3 = simplify(simplify(flux_3));

%Define remaining numerical values.
kon_n = 10^(-3); % diffusion-limited on-rate

% Plug in all numerical values. DB = 0.1*D.
% x = 1/D; y = KD.
flux_final(x,y) = D_n*simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), DB, 0.1*D), ll, ll_n),...
    D, 1/x),koff,kon*y),kon,kon_n));

% flux_final = simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
%     flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
%     DB, 5.88),koff,0.1),kon,kon_n));

% Find the flux when there's no binding. (Should be a constant)
flux_no_binding = D_n*simplify(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    koff, koff),kon,0));

% Normalize the flux to the non-binding flux.
flux_ratio(x,y) = flux_final(x,y)/flux_no_binding;

% Make an x-axis for the plot (do not start at zero; it gets mad).
x_axis = logspace(-1,1);

close all
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 700, 600]);
% Plot flux ratio vs Dtilde/D for several KD values
for kd_n=[0.01, 0.1,1,10] % in uM
flux_ratio_final(x) = flux_ratio(x,kd_n); %#ok<SAGROW>
y_values = double(flux_ratio_final(x_axis));
semilogx(x_axis,y_values,'LineWidth',3);
xlabel('1/D_F (\mus/nm^2)','fontsize',fontsize);
ylabel('Selectivity','fontsize',fontsize);
% ax=gca;
% ax.XTick = [0.1 0.5 1 5 10];
% ax.YTick = [0 50 100 150];
% ax.XTickLabel = {[ fs ' 0.1'],[ fs ' 0.5'],[ fs ' 1'],[ fs ' 5'],[ fs ' 10']};
% ax.YTickLabel = {[ fs ' 0'],[ fs ' 50'],[ fs ' 100'],[ fs ' 150']};
hold all
end
% h = legend({'10 nM','100 nM', '1 uM', '10 uM'});
% title(h,'K_D')
% title('Normalized flux vs. 1/D')

%% Flux ratio vs L for fixed DB, at several KDs
% For Fig. 1 of paper (no functional form of Dtilde)

% Re-dimensionalize without functional form of DB
flux_3 = subs(subs(subs(flux_2, DB, (DB/(koff))),...
    gam, kon*Nt/koff), DF, D/(koff));
flux_3 = simplify(simplify(flux_3));

%Define remaining numerical values.
kon_n = 10^(-3); % diffusion-limited on-rate

% Plug in all numerical values. DB = 0.1*D.
% x = L; y = KD.
flux_final(x,y) = D_n*simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, x), DB, 0.1*D), ll, ll_n),...
    D, D_n),koff,kon*y),kon,kon_n));

% flux_final = simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
%     flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
%     DB, 5.88),koff,0.1),kon,kon_n));

% Find the flux when there's no binding. (Should be a constant)
flux_no_binding(x) = D_n*simplify(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, x), D, D_n), ll, ll_n),...
    koff, koff),kon,0));

% Normalize the flux to the non-binding flux.
flux_ratio(x,y) = flux_final(x,y)/flux_no_binding(x);

% Make an x-axis for the plot (do not start at zero; it gets mad).
x_axis = logspace(1,3);

close all

% Plot flux ratio vs Dtilde/D for several KD values
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 700, 600]);
for kd_n=[0.01, 0.1,1,10] % in uM
flux_ratio_final(x) = flux_ratio(x,kd_n); %#ok<SAGROW>
y_values = double(flux_ratio_final(x_axis));
semilogx(x_axis,y_values,'LineWidth',3);
xlabel('L (nm)','fontsize',fontsize);
ylabel('Selectivity','fontsize',fontsize);
% ax=gca;
% ax.XTick = [10 50 100 500 1000];
% ax.YTick = [100 200 300 400 500];
% ax.XTickLabel = {[ fs ' 10'],[ fs ' 50'],[ fs ' 100'],[ fs ' 500'],[ fs ' 1000']};
% ax.YTickLabel = {[ fs ' 100'],[ fs ' 200'],[ fs ' 300'],[ fs ' 400'],[ fs ' 500']};
hold all
end
% h = legend({'10 nM','100 nM', '1 uM', '10 uM'});
% title(h,'K_D')
% title('Normalized flux vs. L')

%% Absolute flux vs kon for fixed DB, at several KDs
% For Fig. 1 of paper (no functional form of Dtilde)

% Re-dimensionalize without functional form of DB
flux_3 = subs(subs(subs(flux_2, DB, (DB/(koff))),...
    gam, kon*Nt/koff), DF, D/(koff));
flux_3 = simplify(simplify(flux_3));

%Define remaining numerical values.
kon_n = 10^(-3); % diffusion-limited on-rate

% Plug in all numerical values. DB = 0.1*D.
% x = kon; y = KD.
flux_final(x,y) = D_n*simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    DB, 0.1*D_n),koff,kon*y),kon,x));

% flux_final = simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
%     flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
%     DB, 5.88),koff,0.1),kon,kon_n));

% Find the flux when there's no binding. (Should be a constant)
flux_no_binding = D_n*simplify(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    koff, koff),kon,0));

% Normalize the flux to the non-binding flux.
flux_ratio(x,y) = flux_final(x,y)/flux_no_binding;

% Make an x-axis for the plot (do not start at zero; it gets mad).
x_axis = logspace(-5,-3);

close all
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 700, 600]);
figure('DefaultAxesFontSize',18)
% Plot flux ratio vs Dtilde/D for several KD values
for kd_n=[0.01, 0.1,1,10] % in uM
flux_ratio_final(x) = flux_final(x,kd_n); %#ok<SAGROW>
y_values = double(flux_ratio_final(x_axis));
semilogx(x_axis,y_values,'LineWidth',3);
xlabel('k_{on} (\mus^{-1}\muM^{-1})','fontsize',fontsize);
ylabel('Flux (\muM/\mus nm)','fontsize',fontsize);
%  ax=gca;
%  ax.XTick = [1e-5 1e-4 1e-3];
%  ax.YTick = [0.1 0.2 0.3 0.4 0.5];
%  ax.XTickLabel = {[ fs ' 10^{-5}'],[ fs ' 10^{-4}'],[ fs ' 10^{-3}']};
%  ax.YTickLabel = {[ fs ' 0.1'],[ fs ' 0.2'],[ fs ' 0.3'],[ fs ' 0.4'],[ fs ' 0.5']};
hold all
end
% h = legend({'10 nM','100 nM', '1 uM', '10 uM'});
% title(h,'K_D')
% title('Flux vs. k_{on}')

%% Selectivity vs kon for fixed DB, at several KDs
% For Fig. 1 of paper (no functional form of Dtilde)

% Re-dimensionalize without functional form of DB
flux_3 = subs(subs(subs(flux_2, DB, (DB/(koff))),...
    gam, kon*Nt/koff), DF, D/(koff));
flux_3 = simplify(simplify(flux_3));

%Define remaining numerical values.
kon_n = 10^(-3); % diffusion-limited on-rate

% Plug in all numerical values. DB = 0.1*D.
% x = kon; y = KD.
flux_final(x,y) = D_n*simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    DB, 0.1*D_n),koff,kon*y),kon,x));

% flux_final = simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
%     flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
%     DB, 5.88),koff,0.1),kon,kon_n));

% Find the flux when there's no binding. (Should be a constant)
flux_no_binding = D_n*simplify(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    koff, koff),kon,0));

% Normalize the flux to the non-binding flux.
flux_ratio(x,y) = flux_final(x,y)/flux_no_binding;

% Make an x-axis for the plot (do not start at zero; it gets mad).
x_axis = logspace(-5,-3);

close all
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 700, 600]);
figure('DefaultAxesFontSize',18)
% Plot flux ratio vs Dtilde/D for several KD values
for kd_n=[0.01, 0.1,1,10] % in uM
flux_ratio_final(x) = flux_ratio(x,kd_n); 
y_values = double(flux_ratio_final(x_axis));
semilogx(x_axis,y_values,'LineWidth',3);
xlabel('k_{on} (\mus^{-1}\muM^{-1})','fontsize',fontsize);
ylabel('Selectivity','fontsize',fontsize);
% ax=gca;
% ax.XTick = [1e-5 1e-4 1e-3];
% ax.YTick = [10 20 30 40 50];
% ax.XTickLabel = {[ fs ' 10^{-5}'],[ fs ' 10^{-4}'],[ fs ' 10^{-3}']};
% ax.YTickLabel = {[ fs ' 10'],[ fs ' 20'],[ fs ' 30'],[ fs ' 40'],[ fs ' 50']};
hold all
end
% h = legend({'10 nM','100 nM', '1 uM', '10 uM'});
% title(h,'K_D')
% title('Flux vs. k_{on}')

%% Absolute flux vs Nt for fixed DB, at several KDs
% For Fig. 1 of paper (no functional form of Dtilde)

% Re-dimensionalize without functional form of DB
flux_3 = subs(subs(subs(flux_2, DB, (DB/(koff))),...
    gam, kon*Nt/koff), DF, D/(koff));
flux_3 = simplify(simplify(flux_3));

%Define remaining numerical values.
kon_n = 10^(-3); % diffusion-limited on-rate

% Plug in all numerical values. DB = 0.1*D.
% x = Nt; y = KD.
flux_final(x,y) = D_n*simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, x),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    DB, 0.1*D_n),koff,kon*y),kon,kon_n));

% flux_final = simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
%     flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
%     DB, 5.88),koff,0.1),kon,kon_n));

% Find the flux when there's no binding. (Should be a constant)
flux_no_binding = D_n*simplify(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    koff, koff),kon,0));

% Normalize the flux to the non-binding flux.
flux_ratio(x,y) = flux_final(x,y)/flux_no_binding;

% Make an x-axis for the plot (do not start at zero; it gets mad).
x_axis = logspace(1,4);

close all
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 700, 600]);
% Plot flux ratio vs Dtilde/D for several KD values
for kd_n=[0.01, 0.1,1,10] % in uM
flux_ratio_final(x) = flux_final(x,kd_n); %#ok<SAGROW>
y_values = double(flux_ratio_final(x_axis));
semilogx(x_axis,y_values,'LineWidth',3);
xlabel('N_t (\muM)','fontsize',fontsize);
ylabel('Flux (\muM/\mus nm)','fontsize',fontsize);
% ax=gca;
% ax.XTick = [10 100 10^3 10^4];
% ax.YTick = [0.5 1.0 1.5];
% ax.XTickLabel = {[ fs ' 10'],[ fs ' 100'],[ fs ' 10^{3}'],[ fs ' 10^{4}']};
% ax.YTickLabel = {[ fs ' 0.5'],[ fs ' 1.0'],[ fs ' 1.5']};
hold all
end
% h = legend({'10 nM','100 nM', '1 uM', '10 uM'});
% title(h,'K_D')
% title('Flux vs. N_t')

%% Absolute flux vs 1/D for fixed DB, at several KDs
% For Fig. 1 of paper (no functional form of Dtilde)

% Re-dimensionalize without functional form of DB
flux_3 = subs(subs(subs(flux_2, DB, (DB/(koff))),...
    gam, kon*Nt/koff), DF, D/(koff));
flux_3 = simplify(simplify(flux_3));

%Define remaining numerical values.
kon_n = 10^(-3); % diffusion-limited on-rate

% Plug in all numerical values. DB = 0.1*D.
% x = 1/D; y = KD.
flux_final(x,y) = (1/x)*simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), DB, 0.1*D), ll, ll_n),...
    D, 1/x),koff,kon*y),kon,kon_n));

% flux_final = simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
%     flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
%     DB, 5.88),koff,0.1),kon,kon_n));

% Find the flux when there's no binding. (Should be a constant)
flux_no_binding = simplify(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
    koff, koff),kon,0));

% Normalize the flux to the non-binding flux.
flux_ratio(x,y) = flux_final(x,y)/flux_no_binding;

% Make an x-axis for the plot (do not start at zero; it gets mad).
x_axis = logspace(-1,1);

close all
FigHandle = figure('DefaultAxesFontSize',18);
set(FigHandle, 'Position', [100, 100, 700, 600]);
% Plot flux ratio vs Dtilde/D for several KD values
for kd_n=[0.01, 0.1,1,10] % in uM
flux_final2(x) = flux_final(x,kd_n); %#ok<SAGROW>
y_values = double(flux_final2(x_axis));
semilogx(x_axis,y_values,'LineWidth',3);
xlabel('1/D_F (\mus/nm^2)','fontsize',fontsize);
ylabel('Flux (\muM/\mus nm)','fontsize',fontsize);
% ax=gca;
% ax.XTick = [1e-1 1 10];
% %ax.YTick = [1 2 3 4 5];
%  ax.XTickLabel = {[ fs ' 0.1'],[ fs ' 1'],[ fs ' 10']};
% % ax.YTickLabel = {[ fs ' 1'],[ fs ' 2'],[ fs ' 3'],[ fs ' 4'],[ fs ' 5']};
hold all
end
% h = legend({'10 nM','100 nM', '1 uM', '10 uM'});
% title(h,'K_D')
% title('Flux vs. 1/D')

%% Absolute flux vs length for fixed DB, at several KDs
% For Fig. 1 of paper (no functional form of Dtilde)

% Re-dimensionalize without functional form of DB
flux_3 = subs(subs(subs(flux_2, DB, (DB/(koff))),...
    gam, kon*Nt/koff), DF, D/(koff));
flux_3 = simplify(simplify(flux_3));

%Define remaining numerical values.
kon_n = 10^(-3); % diffusion-limited on-rate

% Plug in all numerical values. DB = 0.1*D.
% x = L; y = KD.
flux_final(x,y) = D_n*simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, x), DB, 0.1*D), ll, ll_n),...
    D, D_n),koff,kon*y),kon,kon_n));

% flux_final = simplify(subs(subs(subs(subs(subs(subs(subs(subs(...
%     flux_3, Nt, Nt_n),Abound, Abound_n), L, L_n), D, D_n), ll, ll_n),...
%     DB, 5.88),koff,0.1),kon,kon_n));

% Find the flux when there's no binding. (Should be a constant)
flux_no_binding(x) = D_n*simplify(subs(subs(subs(subs(subs(subs(subs(...
    flux_3, Nt, Nt_n),Abound, Abound_n), L, x), D, D_n), ll, ll_n),...
    koff, koff),kon,0));

% Normalize the flux to the non-binding flux.
flux_ratio(x,y) = flux_final(x,y)/flux_no_binding(x);

% Make an x-axis for the plot (do not start at zero; it gets mad).
x_axis = logspace(1,3);

close all
FigHandle = figure;
set(FigHandle, 'Position', [100, 100, 700, 600]);
% Plot flux ratio vs Dtilde/D for several KD values
for kd_n=[0.01, 0.1,1,10] % in uM
flux_ratio_final(x) = flux_final(x,kd_n); %#ok<SAGROW>
y_values = double(flux_ratio_final(x_axis));
semilogx(x_axis,y_values,'LineWidth',3);
xlabel('L (nm)','fontsize',fontsize);
ylabel('Flux (\muM/\mus nm)','fontsize',fontsize);
% ax=gca;
% ax.XTick = [10 100 1000];
% ax.YTick = [0.1 0.2 0.3 0.4 0.5];
% ax.XTickLabel = {[ fs ' 10'],[ fs ' 100'],[ fs ' 1000']};
% ax.YTickLabel = {[ fs ' 0.1'],[ fs ' 0.2'],[ fs ' 0.3'],[ fs ' 0.4'],[ fs ' 0.5']};
hold all
end
% h = legend({'10 nM','100 nM', '1 uM', '10 uM'});
% title(h,'K_D')
% title('Flux vs. L')

