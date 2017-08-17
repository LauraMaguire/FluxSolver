function [flux_2,varlist] = fluxSolver()

clear all %#ok<CLALL>

% Define symbolic variables (all real numbers)
syms DF DB KD KA Abound Bt gam kon koff koff2 alph ll L lam real
syms b d e f g m x y real
% DF = free diffusion coefficient, non-dim (DF = D/(off))
% DB = bound diffusion coefficient, non-dim (DB = Dtilde/(koff))
% KD = dissociation constant koff/kon, dim.
% KA = association constant kon/koff, dim.
% Abound = Conc. of free NTR at edge of gel, non-dim (Abound = A/Bt)
% Bt = total (free+bound) concentration of Nups, dim.
% gam = kon Bt/koff, non-dim.; useful parameter but not physical
% kon = on-rate, dim.
% koff = off-rate (based on affinity, koff = kon*KD), dim.
% koff2 = off-rate within gel with sliding between FGs, dim.
% alph = ll*koff2, contour length*pers.length*koff2, dim.
% ll = tether contour length times persistence length, dim.
% L = length of pore region, dim.
% b d e f g m x y are constants used at intermediate steps in the solution.

% Make physical assumptions about variables
assume(gam>0)
assume(kon>0)
assume(koff>0)
assume(DF>0)
assume(DB>0)
assume(alph>0)

% Define symbolic functions
syms A(x) B(x) C(x) cx(x)
% A(x) is the concentration of free NTR, non-dim. (A = [A]/Bt)
% B(x) is the concentration of free Nup, non-dim. (B = [B]/Bt)
% C(x) is the concentration of NTR-Nup complex, non-dim. (C = [C]/Bt)
% cx(x) is the portion of C(x) not due to chemical equilibrium, non-dim.


%% Linear Case (Check chemical equilibrium)

% Original binding-diffusion equations (at steady-state, LHS = 0)
eq1_1 = -gam*A + C + DF*diff(A,x,2);
eq2_1 = gam*A - C + DB*diff(C,x,2);

% Make the change of variables C = ce + cx = gam*A + cx
eq1_2 = simplify(subs(eq1_1,C, cx+gam*A));
eq2_2 = simplify(simplify(subs(eq2_1,C, cx+gam*A)));

% Check: Eq. 1 at steady-state should now read:
% 0 = cx + DF*diff(A,x,2)
eq1_3 = simplify(subs(eq1_2, cx, -DF*diff(A(x), x, x)));
% Print an error message if eq1_3 does not equal zero.
if eq1_3 ~= 0
    fprintf('I am broken at eq1_3.') 
end 

% Substitute in the expression for cx found using eq1_3 above:
% cx = -DF*diff(A,x,2)
eq2_3 = simplify(simplify(subs(eq2_2, cx, -DF*diff(A(x), x, x))));
% Next we need to solve eq2_3 = 0 (because of steady-state).

%% Linear case solutions
% Linear case uses the assumption that most Nups are not bound (Bt>>C)

% Define the dimensionless parameter lam:
%lam = sqrt((DF+DB*gam)/(DF*DB));
% Guess a solution for A(x) in the linear case:
% (Parameters m, b, f, and g must be deterimined using BCs.)
A_1(x) = m*x + b+ f*exp(x*lam)+g*exp(-lam*x);

% Check that the guess is consistent with the DE:
Check_1 = simplify(subs(eq2_3, A(x), A_1));
% Print an error message if Check_1 does not equal zero.
if Check_1 ~= 0
    if simplify(subs(Check_1, lam, sqrt((DF+DB*gam)/(DF*DB)))) ~=0
    fprintf('I am broken at Check_1.') 
    end
end

% Back-substitute to find other functions in terms of m, b, f, and g
cx_1(x) = simplify(-DF*diff(A_1(x), x, x));
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
f_temp = simplify(solve(subs(difC_1(L), g, g_temp),f));

% The right hand boundary gives m.
% Substitute 'temp' parameters into A_1
A_2(x) = simplify(subs(subs(subs(A_1(x),b,b_temp), g, g_temp),f,f_temp));
% Solve for m:
m_temp = simplify(simplify(solve(A_2(L)-Abound, m)));


%% Finish substituting into the guess solutions to remove unknown parameters.

% Only m needs to be removed from A_temp.
A_final(x) = simplify(subs(A_2(x), m, m_temp));

% Use final solution A_2 to back-substitute and find other solutions.
cx_final(x) = simplify(-DF*diff(A_final(x), x, x));
C_final(x) = cx_final + gam*A_final;
difC_final(x) = simplify(diff(C_final,x));

%% Checks
%Check A+C is linear (only works if DB=DF)
Check_2 = simplify(A_final + C_final);
%Check_2_sub = simplify(subs(subs(Check_2, DB, DF), lam,sqrt((DF+DF*gam)/(DF*DF))));
Check_2_sub = simplify(subs(subs(Check_2, lam, sqrt((DF+DB*gam)/(DF*DB))), DB, DF));

% Check that the flux of C at each edge is zero.
Check_3a = simplify(difC_final(0));
if Check_3a ~= 0
    fprintf('I am broken at Check_3a.') 
end

Check_3b = simplify(difC_final(L));
if Check_3b ~= 0
    fprintf('I am broken at Check_3b.') 
end

% Check that the concentration of A at each edge is correct.
Check_4a = A_final(0);
if Check_4a ~= 0
    fprintf('I am broken at Check_4a.') 
end

Check_4b = simplify(A_final(L));
if (Check_4b - Abound) ~= 0
    fprintf('I am broken at Check_4b.') 
end


%% Determine the flux out of the pore using Fick's law.

% Take the spatial derivative of A_final(x).
flux_1 = simplify(subs(diff(A_final,x),lam,sqrt((DF+DB*gam)/(DF*DB))));

% Set x=0 to find the flux at the edge of the pore.
flux_2 = simplify(flux_1(0));

%% Set reasonable parameter values.

% alph is lc*lp*koff =ll*koff
% gam = kon*Nt/koff 
% lam = sqrt((DF+DB*gam)/(DF*DB));
% DF is D/(L^2*koff) 
% Dtilde  = bound state diffusion =(D*alph)/(alph+3*D); 
% alph is lc*lp*koff (propeties of the polymer)
% DB = Dtilde/L^2*koff

% Define Nt (same as Bt) and D (free diffusion coefficient with dimensions)
syms Nt D

% Set variables to reasonable values (SI units but in us, nm, and uM).
Nt_n = 10^3; %total nup concentration in uM
Abound_n = 1; % 1 % free NTR at left boundary in uM (upper bound)
L_n = 100; % nm, length of pore
D_n = 1; %nm^/us INSIDE THE PORE
ll_n = 250; % nm^2 (1 nm)*(countour length, 1 ammino acid is roughly 1 nm).

% Clear and define dummy variables for plotting.
clear x y
syms x y


%% Start filling in actual numbers in place of variables.

fluxFinal = simplify(subs(flux_2,DF,D_n));
varlist = struct();
varlist.DF = DF;
%syms DF DB KD KA Abound Bt gam kon koff koff2 alph ll L lam real

end