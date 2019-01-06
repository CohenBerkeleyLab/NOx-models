function [oh, ho2, ro2, vocr_soln, tau_soln] = hox_solve_tau_constraint(nox, phox, tau, alpha, varargin)
% HOX_SOLVE_TAU_CONSTRAINT Solve for OH, HO2, RO2, and VOCR given NOx, PHOx, alpha, and total NOx lifetime.
%
%   [ OH, HO2, RO2, VOCR ] = HOX_SOLVE_TAU_CONSTRAINT( NOX, PHOX, TAU,
%   ALPHA ) Given [NOx] (in molec. cm^-3), P(HOx) (production of HOx in
%   molec. cm^-3 s^-1), total NOx lifetime (tau, in hours), and the RO2 +
%   NO branching ratio (alpha, unitless), attempts to solve for the steady
%   state concentrations of OH, HO2, and RO2 as well as the VOC reactivity
%   (VOCR) that are consistent with the given constraints.
%
%   Parameters:
%
%       'fmincon_debug' - controls how much printing to the console fmincon
%       does. Default is 'none', other options are 'final' and 'iter'.

p = advInputParser;
p.addParameter('fmincon_debug', 'none');
p.parse(varargin{:});
pout = p.Results;

fmincon_debug = pout.fmincon_debug;

cost_fxn = @(vocr) (nox_lifetime(nox, 'phox', phox, 'alpha', alpha, 'vocr', vocr) - tau).^2;

history.x = [];
opts = optimoptions('fmincon','Display',fmincon_debug,'OutputFcn',@outfun);

vocr_soln = fmincon(cost_fxn, 1, -1, 0,[],[],[],[],[], opts);%, -1, 0);
[oh, ho2, ro2] = hox_ss_solver(nox, phox, vocr_soln, alpha);
tau_soln = nox_lifetime(nox, 'phox', phox, 'alpha', alpha, 'vocr', vocr_soln);

    function stop = outfun(x,optimvals,state)
        stop = false;
        switch state
            case 'iter'
                history.x = [history.x; x];
        end
    end
end

