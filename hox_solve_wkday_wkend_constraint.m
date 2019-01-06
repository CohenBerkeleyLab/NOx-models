function [wkday, wkend, history] = hox_solve_wkday_wkend_constraint(nox_ratio, tau_wkday, tau_wkend, phox, alpha, varargin)
% HOX_SOLVE_WKDAY_WKEND_CONSTRAINT Solve for OH, HO2, RO2, and VOCR from a weekend/weekday NOx ratio and lifetimes.
%
%   [ OH, HO2, RO2, VOCR ] = HOX_SOLVE_WKDAY_WKEND_CONSTRAINT( NOX_RATIO,
%   TAU_WKDAY, TAU_WKEND, PHOX, ALPHA)
%
%       Inputs are:
%           * NOX_RATIO - the ratio of weekend NOx/weekday NOx
%           * TAU_WKDAY, TAU_WKEND - the weekday and weekend NOx lifetimes
%             in hours
%           * PHOX - the production of HOx in molec. cm^-3 s^-1
%           * ALPHA - the RO2 + NO branching ratio (unitless)
%
%       This will attempt to solve for the weekend and weekday NOx, OH,
%       HO2, and RO2 concentrations, as well as VOCR value by searching for
%       a pair of points on the NOx lifetime vs. [NOx] and VOCR isopleth
%       that satisfies the given weekend/weekday NOx ratio and the given
%       weekend and weekday lifetimes. PHOX and ALPHA are used in the HOx
%       steady state solution.
%
%       Returned structures have nox, oh, ho2, ro2 concentrations in molec.
%       cm^-3 and vocr in s^-1.
%
%   Parameters:
%       'nox_initial' - initial NOx concentration, in ppb, to use in the
%       solver.
%
%       'fmincon_debug' - controls how much printing to the console fmincon
%       does. Default is 'none', other options are 'final' and 'iter'.

p = advInputParser;
p.addParameter('nox_initial', 5);
p.addParameter('fmincon_debug', 'none');
p.parse(varargin{:});
pout = p.Results;

nox_initial = pout.nox_initial;
fmincon_debug = pout.fmincon_debug;

% It seems like fmincon struggles when the two variables have very
% different orders of magnitude - when [NOx] is ~10^11 and VOCR is ~10, it
% seems to only want to vary VOCR by any appreciable percentage. The
% solution was to let fmincon vary [NOx] in ppb, but in the cost function
% transform it to molec. cm^-3 for the lifetime calculations, and tranform
% the final result back to number density, since that what was actually
% used in the minimization.
ppb2numdens = 2e10;

history.x = [];
opts = optimoptions('fmincon','Display',fmincon_debug,'OutputFcn',@outfun);

% assume 5 ppb weekday NOx and VOCR = 1 to start
x0 = [nox_initial, 1]; 
% Ax <= b: require that NOx >= 0 and VOCR >= 0.5. 
A = diag([-1 -1]); 
b = [0; 0];
%b = [0; -0.5];
warning('off', 'symbolic:numeric:NumericalInstability')
soln = fmincon(@cost_function, x0, A, b,[],[],[],[],[], opts);%, -1, 0);
warning('on', 'symbolic:numeric:NumericalInstability');
nox_soln = soln(1) * ppb2numdens;
vocr_soln = soln(2);

[oh, ho2, ro2] = hox_ss_solver(nox_soln, phox, vocr_soln, alpha);
tau_soln = nox_lifetime(nox_soln, 'phox', phox, 'alpha', alpha, 'vocr', vocr_soln);
wkday = struct('nox', nox_soln, 'oh', oh, 'ho2', ho2, 'ro2', ro2, 'vocr', vocr_soln, 'tau', tau_soln);

[oh, ho2, ro2] = hox_ss_solver(nox_soln * nox_ratio, phox, vocr_soln, alpha);
tau_soln = nox_lifetime(nox_soln * nox_ratio, 'phox', phox, 'alpha', alpha, 'vocr', vocr_soln);
wkend = struct('nox', nox_soln * nox_ratio, 'oh', oh, 'ho2', ho2, 'ro2', ro2, 'vocr', vocr_soln, 'tau', tau_soln);


    function stop = outfun(x,optimvals,state)
        stop = false;
        switch state
            case 'iter'
                history.x = [history.x; x];
        end
    end

    function cost = cost_function(x)
        wkday_nox = x(1) * ppb2numdens;
        vocr = x(2);
        model_tau_wkday = nox_lifetime(wkday_nox, 'phox', phox, 'alpha', alpha, 'vocr', vocr);
        model_tau_wkend = nox_lifetime(wkday_nox * nox_ratio, 'phox', phox, 'alpha', alpha, 'vocr', vocr);
        cost = (tau_wkday - model_tau_wkday).^2 + (tau_wkend - model_tau_wkend).^2;
    end
end

