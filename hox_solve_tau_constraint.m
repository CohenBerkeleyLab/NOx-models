function [oh, ho2, ro2, soln] = hox_solve_tau_constraint(nox, phox, tau, alpha, varargin)
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
p.addParameter('variables', {'vocr'});
p.addParameter('fmincon_debug', 'none');
p.parse(varargin{:});
pout = p.Results;

E = JLLErrors;

variables = pout.variables;
fmincon_debug = pout.fmincon_debug;

vocr_init = 1;
phox_convert = 2e19*1e-12;
n_vars = numel(variables);
lb = zeros(n_vars,1);
ub = Inf(n_vars,1);

% this is a kludge to make sure that the variables are in the right order
% for when we put the results in the output struct
variables = fliplr(sort(variables));

if all(ismember({'vocr','phox','alpha'}, variables))
    cost_fxn = @cost_fxn_all;
    x0 = [vocr_init; phox ./ phox_convert; alpha];
    ub(3) = 1;
elseif all(ismember({'vocr','phox'}, variables))
    cost_fxn = @cost_fxn_vocr_phox;
    x0 = [vocr_init; phox ./ phox_convert];
elseif all(ismember({'vocr','alpha'}, variables))
    cost_fxn = @cost_fxn_vocr_alpha;
    x0 = [vocr_init; alpha];
    ub(2) = 1;
elseif all(ismember({'vocr'}, variables))
    cost_fxn = @cost_fxn_vocr;
    x0 = [vocr_init];
else
    E.notimplemented('No cost function implemented for the combination of variables: %s', strjoin(variables, ', '));
end

%A = diag(-ones(n_vars,1));
%b = zeros(n_vars,1);


history.x = [];
opts = optimoptions('fmincon','Display',fmincon_debug,'OutputFcn',@outfun);

fmin_out = fmincon(cost_fxn, x0, [], [], [], [], lb, ub, [], opts);%, -1, 0);
soln = struct('vocr', [], 'phox', phox ./ phox_convert, 'alpha', alpha, 'tau', []);
for i = 1:n_vars
    soln.(variables{i}) = fmin_out(i);
end
soln.phox = soln.phox .* phox_convert;

[oh, ho2, ro2] = hox_ss_solver(nox, soln.phox, soln.vocr, soln.alpha);
soln.tau = nox_lifetime(nox, 'phox', soln.phox, 'alpha', soln.alpha, 'vocr', soln.vocr);

    function stop = outfun(x,optimvals,state)
        stop = false;
        switch state
            case 'iter'
                history.x = [history.x; x];
        end
    end

    function cost = cost_fxn_vocr(vocr)
        cost = nox_lifetime(nox, 'phox', phox, 'alpha', alpha, 'vocr', vocr) - tau;
        cost = cost .^ 2;
    end

    function cost = cost_fxn_vocr_phox(x)
        vocr = x(1);
        phox = x(2) .* phox_convert;
        cost = nox_lifetime(nox, 'phox', phox, 'alpha', alpha, 'vocr', vocr) - tau;
        cost = cost .^ 2;
    end

    function cost = cost_fxn_vocr_alpha(x)
        vocr = x(1);
        alpha = x(2);
        cost = nox_lifetime(nox, 'phox', phox, 'alpha', alpha, 'vocr', vocr) - tau;
        cost = cost .^ 2;
    end

    function cost = cost_fxn_all(x)
        vocr = x(1);
        phox = x(2) .* phox_convert;
        alpha = x(3);
        cost = nox_lifetime(nox, 'phox', phox, 'alpha', alpha, 'vocr', vocr) - tau;
        cost = cost .^ 2;
    end
end

