function [bindFlux, nonbindFlux] = subNum(params,tetherFlag)
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
syms symflux gam alph real
syms DF AB Nt L DBb kon koff ll real

% Begin with fully-symbolic expression for flux.
%symflux = (AB*(DF+DBb*gam)^(3/2)*(exp((L*(DF+DBb*gam)^(1/2))/(DBb^(1/2)*DF^(1/2)))+1))/...
%    (L*(DF + DBb*gam)^(3/2) - 2*DBb^(3/2)*DF^(1/2)*gam + ...
%    L*exp((L*(DF+DBb*gam)^(1/2))/(DBb^(1/2)*DF^(1/2)))*(DF+DBb*gam)^(3/2)-...
%    DBb*L*gam*(DF + DBb*gam)^(1/2) + ...
%    2*DBb^(3/2)*DF^(1/2)*gam*exp((L*(DF+DBb*gam)^(1/2))/(DBb^(1/2)*DF^(1/2)))-...
%    DBb*L*gam*exp((L*(DF + DBb*gam)^(1/2))/(DBb^(1/2)*DF^(1/2)))*(DF + DBb*gam)^(1/2));

symflux = (AB*(DF + DBb*gam)^(3/2)*(exp((L*(DF + DBb*gam)^(1/2))/(DBb^(1/2)*DF^(1/2))) + 1))/(L*(DF + DBb*gam)^(3/2) - 2*DBb^(3/2)*DF^(1/2)*gam + L*exp((L*(DF + DBb*gam)^(1/2))/(DBb^(1/2)*DF^(1/2)))*(DF + DBb*gam)^(3/2) - DBb*L*gam*(DF + DBb*gam)^(1/2) + 2*DBb^(3/2)*DF^(1/2)*gam*exp((L*(DF + DBb*gam)^(1/2))/(DBb^(1/2)*DF^(1/2))) - DBb*L*gam*exp((L*(DF + DBb*gam)^(1/2))/(DBb^(1/2)*DF^(1/2)))*(DF + DBb*gam)^(1/2));

% Replace gam with kon*Nt/koff
symflux = subs(symflux, gam, kon*Nt/koff);

% Put flux expression fully into dimensional terms
symflux = subs(symflux, DBb, DBb/koff);
symflux = simplify(subs(symflux, DF, DF/koff));

% Substitute numerical values for AB, L, and Nt.
symflux = subs(symflux, AB, params.AB);
symflux = subs(symflux, L, params.L);
symflux = simplify(subs(symflux, Nt, params.Nt));

% tetherFlag determines whether tethered-model expression for DBb is used.
if tetherFlag ==1
    % Replace DBb with analytic expression from tether model.
    symflux = subs(symflux, DBb, (DF*alph)/(alph+3*DF));
    symflux = subs(symflux, alph, ll*koff);
    % Substitute a numerical value for ll
    symflux = simplify(subs(symflux, ll, params.ll));
    
else
    % Directly substitute a numerical value for DBb
    symflux = simplify(subs(symflux, DBb, params.DB));
end

% Finally, replace koff and DF.
symflux = subs(symflux, koff, params.koff);
symflux = subs(symflux, DF, params.DF);

% Flux with binding: replace kon with numerical value.
symfluxb = simplify(subs(symflux, kon, params.kon)); % substitute for kon

% Flux without binding: set kon = 0.
symfluxn = simplify(subs(symflux, kon, 0)); % set kon = 0 (no binding)

% Return the numerical flux.
bindFlux = double(symfluxb);
nonbindFlux = double(symfluxn);

end