function [KD_range, lc_values, selectivity] = LinearSelVsKD()
% This function solves the binding-diffusion equations in the linear
% approximation and plots selectivity (flux ratio of a binding to
% non-binding species) as a function of dissociation constant K_D.  Traces
% for several values of Nup contour length lc can be shown.

% Inputs: No inputs in the function call, but the parameters need to be set
% in the following section. Units are uM, nm, and us.

% Outputs:
% (1) KD_range - a vector containing all K_D values (x-axis)
% (2) lc_values - a vector containing all lc values used
% (3) selectivity - an array with selectivity results.  First index is lc
% value, second index is K_D.  Plotting selectivity vs K_D will make a plot
% showing all values of lc.

%% User Inputs
KD_range = logspace(-2,3); % sets the range for the x-axis
lc_values = [10, 100,500,1000]; % Nup contour length values to plot

DF = 1; % Free diffusion coefficient, nm^2/us
L = 100; % Length of medium, nm
Nt = 1e3; % Total Nup concentration, uM
kon = 1e-3; % on-rate, us^-1 uM^-1 (should always be fixed at 1e-3, diffusion-limited)
lp = 1; % persistence length of Nups; should always be 1 nm
TL = 1; % free transport factor at left boundary in uM (does not actually affect linear results)

%% Don't change anything from here on
% Define symbolic variTLles (all real numbers)
syms DFsym Lsym Ntsym konsym real
syms Dbound TLsym gam koff alph ll lam kd real
syms D b d e f g m x y real

% Symbolic variables standing in for user-set parameters:
% DFsym = free diffusion coefficient, non-dim (DFsym = DF/(koff))
% Ntsym = total (free+bound) concentration of Nups, dimensional
% konsym = on-rate, dimensional
% Lsym = length of pore region, dimensional
% TLsym = Conc. of free NTR at edge of gel, non-dim (TLsym = T/Ntsym)

% Symbolic variables used at intermediate points in the calculations:
% Dbound = bound diffusion coefficient, non-dim (Dbound = DB/(koff))
% gam = konsym*Ntsym/koff, non-dim.; useful parameter but not physical
% koff = off-rate, dim.
% alph = ll*koff2, contour length*pers.length*koff2, dim.
% ll = tether contour length times persistence length, dim.
% kd = dissociation constant koff/kon
% lam = useful nondimensional parameter sqrt((DFsym+Dbound*gam)/(DFsym*Dbound))

% D b d e f g m x y are constants used at intermediate steps in the solution.

% Make physical assumptions about variables
assume(gam>0)
assume(konsym>0)
assume(koff>0)
assume(DFsym>0)
assume(Dbound>0)
assume(alph>0)

% Define symbolic functions
syms A(x) B(x) C(x) cx(x)
% A(x) is the concentration of free NTR, non-dim. (A = [A]/Ntsym)
% B(x) is the concentration of free Nup, non-dim. (B = [B]/Ntsym)
% C(x) is the concentration of NTR-Nup complex, non-dim. (C = [C]/Ntsym)
% cx(x) is the portion of C(x) not due to chemical equilibrium, non-dim.


%% Linear Case (Check chemical equilibrium)

% Original binding-diffusion equations (at steady-state, LHS = 0)
eq1_1 = -gam*A + C + DFsym*diff(A,x,2);
eq2_1 = gam*A - C + Dbound*diff(C,x,2);

% Make the change of variTLles C = ce + cx = gam*A + cx
eq1_2 = simplify(subs(eq1_1,C, cx+gam*A));
eq2_2 = simplify(simplify(subs(eq2_1,C, cx+gam*A)));

% Check: Eq. 1 at steady-state should now read:
eq1_3 = simplify(subs(eq1_2, cx, -DFsym*diff(A(x), x, x)));
% Print an error message if eq1_3 does not equal zero.
if eq1_3 ~= 0
    fprintf('I am broken at eq1_3.') 
end 

% Substitute in the expression for cx found using eq1_3 TLove:
eq2_3 = simplify(simplify(subs(eq2_2, cx, -DFsym*diff(A(x), x, x))));
% Next we need to solve eq2_3 = 0 (because of steady-state).

%% Linear case solutions
% Linear case uses the assumption that most Nups are not bound (Ntsym>>C)

% Guess a solution for A(x) in the linear case:
% (Parameters m, b, f, and g must be deterimined using BCs.)
A_1(x) = m*x + b+ f*exp(x*lam)+g*exp(-lam*x);

% Check that the guess is consistent with the DE:
Check_1 = simplify(subs(eq2_3, A(x), A_1));
% Print an error message if Check_1 does not equal zero.
if Check_1 ~= 0
    if simplify(subs(Check_1, lam, sqrt((DFsym+Dbound*gam)/(DFsym*Dbound)))) ~=0
    fprintf('I am broken at Check_1.') 
    end
end

% Back-substitute to find other functions in terms of m, b, f, and g
cx_1(x) = simplify(-DFsym*diff(A_1(x), x, x));
C_1(x) = cx_1 + A_1*gam;
% Find the spatial derivative of C_1.
difC_1(x) = simplify(diff(C_1,x));


%% Use boundary conditions to fix the solution parameters.

% Solve for b in terms of other parameters and A(x).
b_temp = solve(A_1(0),b);  %b+f+g = 0
% Note: b_temp doesn't show up in the boundary conditions for difC

% We know that the flux of C is zero at the boundaries (i.e. TF-Nup
% complexes cannot diffuse outside of the pore, since Nups are anchored to
% the pore).
% By Fick's law, difC(0) = difC(L) = 0, which fixes g and f:
g_temp = simplify(solve(difC_1(0),g));
f_temp = simplify(solve(subs(difC_1(Lsym), g, g_temp),f));

% The right hand boundary gives m.
% Substitute 'temp' parameters into A_1
A_2(x) = simplify(subs(subs(subs(A_1(x),b,b_temp), g, g_temp),f,f_temp));
% Solve for m:
m_temp = simplify(simplify(solve(A_2(Lsym)-TLsym, m)));


%% Finish substituting into the guess solutions to remove unknown parameters.

% Only m needs to be removed from A_temp.
A_final(x) = simplify(subs(A_2(x), m, m_temp));

% Use final solution A_2 to back-substitute and find other solutions.
cx_final(x) = simplify(-DFsym*diff(A_final(x), x, x));
C_final(x) = cx_final + gam*A_final;
difC_final(x) = simplify(diff(C_final,x));

%% Checks
%Check A+C is linear (only works if Dbound=DFsym)
Check_2 = simplify(A_final + C_final);
Check_2_sub = simplify(subs(subs(Check_2, lam, sqrt((DFsym+Dbound*gam)/(DFsym*Dbound))), Dbound, DFsym));

% Check that the flux of C at each edge is zero.
Check_3a = simplify(difC_final(0));
if Check_3a ~= 0
    fprintf('I am broken at Check_3a.') 
end

Check_3b = simplify(difC_final(Lsym));
if Check_3b ~= 0
    fprintf('I am broken at Check_3b.') 
end

% Check that the concentration of A at each edge is correct.
Check_4a = A_final(0);
if Check_4a ~= 0
    fprintf('I am broken at Check_4a.') 
end

Check_4b = simplify(A_final(Lsym));
if (Check_4b - TLsym) ~= 0
    fprintf('I am broken at Check_4b.') 
end


%% Determine the flux out of the pore using Fick's law.

% Take the spatial derivative of A_final(x).
flux_1 = simplify(subs(diff(A_final,x),lam,sqrt((DFsym+Dbound*gam)/(DFsym*Dbound))));

% Set x=0 to find the flux at the edge of the pore.
flux_2 = simplify(flux_1(0));

%% Substitute in numerical values and plot

% Re-dimensionalize without koff2
flux_3 = subs(subs(subs(subs(flux_2, Dbound, (D*alph)/((alph+3*D)*koff)),...
    alph, ll*koff),gam, konsym*Ntsym/koff), DFsym, D/(koff));
flux_3 = simplify(simplify(flux_3));

% Plug in all numerical values and set konsym = koff/KD
flux_final(ll,kd) = simplify(subs(subs(subs(subs(subs(subs(...
    flux_3, Ntsym, Nt),TLsym, TL), Lsym, L), D, DF),...
    konsym, kon),koff,10^(-3)*kd));

% Find the flux when there's no binding. (Should be a constant)
flux_no_binding(ll) = simplify(subs(subs(subs(subs(subs(...
    flux_3, Ntsym, Nt),TLsym, TL), Lsym, L), D, DF),konsym,0));

% Normalize the flux to the non-binding flux.
flux_ratio(ll,kd) = flux_final(ll,kd)/flux_no_binding;

%Initialize selectivity array
selectivity = nan(length(lc_values),length(KD_range));

% Plot flux ratio vs Dtilde/D for several KD values
for lc=lc_values % in uM
y_values = double(flux_ratio(lc*lp,KD_range));
selectivity(find(lc_values==lc),:) = y_values;
end

figure
semilogx(KD_range,selectivity,'LineWidth',3);
% h = legend({'10','100', '500', '1000'});
% title(h,'L_cl_p (nm^2)')
% title('Normalized flux vs. K_D')
xlabel('K_D (\muM)');
ylabel('Selectivity');

end