function [numericalA, numericalC] = ACSubNum(params, tetherFlag, bcFlag)
% subNum substitutes numbers into the symbolic flux expression found with
% fluxSolver.
% Inputs:
%   1) params: numbers corresponding to the necessary parameters
%   2) tetherFlag: 1 means "use tethered diffusion expression," 0 means
%   "directly plug in a given value for bound diffusion coefficient."
% Outputs:
%   1) bindFlux: the numerical flux value when binding is allowed.
%   2) nonbindFlux: the numerical flux value when binding is not allowed.

% Define necessary symbolic variables.
syms A C gam lam alph real
syms DF AB Nt L DBb kon koff ll x real

if bcFlag == 0
    A = -(AB*(exp(L*lam) + 1)*(gam*exp(lam*(L + x)) - gam*exp(L*lam) - gam*exp(lam*x) + gam*exp(2*lam*x) + DF*lam^3*x*exp(lam*(L + x)) - gam*lam*x*exp(lam*x) + DF*lam^3*x*exp(lam*x) - gam*lam*x*exp(lam*(L + x))))/((exp(lam*(L + x)) + exp(lam*x))*(2*gam - 2*gam*exp(L*lam) + L*gam*lam - DF*L*lam^3 - DF*L*lam^3*exp(L*lam) + L*gam*lam*exp(L*lam)));
    C = (AB*DF*gam*lam^2*(exp(L*lam) - exp(2*lam*x)))/(2*gam*exp(lam*(L + x)) - 2*gam*exp(lam*x) + DF*L*lam^3*exp(lam*(L + x)) - L*gam*lam*exp(lam*x) + DF*L*lam^3*exp(lam*x) - L*gam*lam*exp(lam*(L + x))) - (AB*gam*(exp(L*lam) + 1)*(gam*exp(lam*(L + x)) - gam*exp(L*lam) - gam*exp(lam*x) + gam*exp(2*lam*x) + DF*lam^3*x*exp(lam*(L + x)) - gam*lam*x*exp(lam*x) + DF*lam^3*x*exp(lam*x) - gam*lam*x*exp(lam*(L + x))))/((exp(lam*(L + x)) + exp(lam*x))*(2*gam - 2*gam*exp(L*lam) + L*gam*lam - DF*L*lam^3 - DF*L*lam^3*exp(L*lam) + L*gam*lam*exp(L*lam)));
else
    A = AB;
    C = gam*AB;
end
    
% Replace lam with sqrt((DF+DB*gam)/(DF*DB))
A = subs(A, lam, sqrt((DF+DBb*gam)/(DF*DBb)));
C = subs(C, lam, sqrt((DF+DBb*gam)/(DF*DBb)));

% Replace gam with kon*Nt/koff
A = subs(A, gam, kon*Nt/koff);
C = subs(C, gam, kon*Nt/koff);

% Put flux expression fully into dimensional terms
A = subs(A, DBb, DBb/koff);
C = subs(C, DBb, DBb/koff);
A = simplify(subs(A, DF, DF/koff));
C = simplify(subs(C, DF, DF/koff));

% Substitute numerical values for AB, L, and Nt.
A = subs(A, AB, params.AB);
A = subs(A, L, params.L);
A = simplify(subs(A, Nt, params.Nt));

C = subs(C, AB, params.AB);
C = subs(C, L, params.L);
C = simplify(subs(C, Nt, params.Nt));

% tetherFlag determines whether tethered-model expression for DBb is used.
if tetherFlag ==1
    % Replace DBb with analytic expression from tether model.
    A = subs(A, DBb, (DF*alph)/(alph+3*DF));
    A = subs(A, alph, ll*koff);
    % Substitute a numerical value for ll
    A = simplify(subs(A, ll, params.ll));
    
    % Replace DBb with analytic expression from tether model.
    C = subs(C, DBb, (DF*alph)/(alph+3*DF));
    C = subs(C, alph, ll*koff);
    % Substitute a numerical value for ll
    C = simplify(subs(C, ll, params.ll));
    
else
    % Directly substitute a numerical value for DBb
    A = simplify(subs(A, DBb, params.DB));
    C = simplify(subs(C, DBb, params.DB));
end

% Replace koff and DF.
A = subs(A, koff, params.koff);
A = subs(A, DF, params.DF);

C = subs(C, koff, params.koff);
C = subs(C, DF, params.DF);

% Replace kon with numerical value.
A = simplify(subs(A, kon, params.kon)); % substitute for kon
C = simplify(subs(C, kon, params.kon)); % substitute for kon

% numericalA = zeros(1,length(params.x));
% numericalC = zeros(1,length(params.x));

    A = simplify(subs(A, x, params.x));
    C = simplify(subs(C, x, params.x));

    numericalA = double(A);
    numericalC = double(C);

% Replace x with numerical value.
% for i = 1:length(params.x)
%     disp(num2str(params.x));
%     A = simplify(subs(A, x, params.x(i)));
%     C = simplify(subs(C, x, params.x(i)));
% 
%     numericalA(i) = double(A);
%     numericalC(i) = double(C);
% end
end

