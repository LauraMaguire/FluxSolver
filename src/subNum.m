function [bindFlux, nonbindFlux] = subNum(params,tetherFlag)

syms symflux gam lam alph real
syms DF AB Nt L DB kon koff ll real

% Begin with fully-symbolic expression for flux.
symflux = (AB*(DF+DB*gam)^(3/2)*(exp((L*(DF+DB*gam)^(1/2))/(DB^(1/2)*DF^(1/2)))+1))/...
    (L*(DF + DB*gam)^(3/2) - 2*DB^(3/2)*DF^(1/2)*gam + ...
    L*exp((L*(DF+DB*gam)^(1/2))/(DB^(1/2)*DF^(1/2)))*(DF+DB*gam)^(3/2)-...
    DB*L*gam*(DF + DB*gam)^(1/2) + ...
    2*DB^(3/2)*DF^(1/2)*gam*exp((L*(DF+DB*gam)^(1/2))/(DB^(1/2)*DF^(1/2)))-...
    DB*L*gam*exp((L*(DF + DB*gam)^(1/2))/(DB^(1/2)*DF^(1/2)))*(DF + DB*gam)^(1/2));

% Replace gam with kon*Nt/koff
symflux = simplify(subs(symflux, gam, kon*Nt/koff));

% Put flux expression fully into dimensional terms
symflux = simplify(subs(symflux, DB, DB/koff));
symflux = simplify(subs(symflux, DF, DF/koff));

% Substitute numerical values for AB, L, and Nt.
symflux = simplify(subs(symflux, AB, params.AB));
symflux = simplify(subs(symflux, L, params.L));
symflux = simplify(subs(symflux, Nt, params.Nt));

% tetherFlag determines whether tethered-model expression for DB is used.
if tetherFlag ==1
    % Replace DB with analytic expression from tether model.
    symflux = simplify(subs(symflux, DB, (DF*alph)/((alph+3*DF)*koff)));
    symflux = simplify(subs(symflux, alph, ll*koff));
    % Substitute a numerical value for ll
    symflux = simplify(subs(symflux, ll, params.ll));
    
else
    % Directly substitute a numerical value for DB
    symflux = simplify(subs(symflux, DB, params.DB));
end

% Finally, replace koff and DF.
symflux = simplify(subs(symflux, koff, params.koff));
symflux = simplify(subs(symflux, DF, params.DF));

% Flux with binding: replace kon with numerical value.
symfluxb = simplify(subs(symflux, kon, params.kon)); % substitute for kon

% Flux without binding: set kon = 0.
symfluxn = simplify(subs(symflux, kon, 0)); % set kon = 0 (no binding)

% Return the numerical flux.
bindFlux = double(symfluxb);
nonbindFlux = double(symfluxn);

end