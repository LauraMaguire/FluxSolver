%% recreate Fig. 1 from the Wyman 66 paper showing dissolved oxygen pressure as a function of position
% For reviewer response to PRE submission of modeling paper, 8/28/19

% Follows Wyman's method of generating the plot (recreates data in Table 1)

p = [2.5 5 10 20 40 100 500 1000]';
y = (p./10).^3./(1+(p./10).^3);
Fx = (2e-12.*y +1e-14.*p);
x_1 = Fx/1.5e-10;
x_3 = Fx/3e-10;
x_6 = Fx/6e-10;

plot(1e2*x_1,p)
hold on
plot(1e2*x_3,p)
plot(1e2*x_6,p)
hold off
xlim([0 2.5])