function [ss_oh,ss_vocr,ss_phox,ss_alpha,ss_tau] = match_wrf_oh(wrf_oh, wrf_nox, tau, varargin)
% MATCH_WRF_OH Solve for the VOCR, PHOx, and ALPHA necessary to match WRF [OH]
%   [SS_OH, SS_VOCR, SS_PHOX, SS_ALPHA, SS_TAU] = MATCH_WRF_OH( WRF_OH,
%   WRF_NOX, TAU) Given the concentration of OH and NOx from WRF (in molec.
%   cm^{-3}) and the lifetime observed from space (in hours), this function
%   calculates what VOCR, PHOx, and ALPHA are necessary to make the steady
%   state model give the same OH concentration for the specified NOx and
%   tau.
%
%   Parameters:
%
%       'show_vals' - prints the values of VOCR, PHOx, and ALPHA at each
%       step of the fmincon solver (true/false).
%
%       'fmincon_debug' - set how much detail fmincon prints
%       ('none'/'final'/'iter')

p = advInputParser;
p.addParameter('show_vals', false);
p.addParameter('fmincon_debug', 'none');
p.parse(varargin{:});
pout = p.Results;

show_vals = pout.show_vals;
fmincon_debug = pout.fmincon_debug;

history.x = [];
opts = optimoptions('fmincon','Display',fmincon_debug,'OutputFcn',@outfun);

% We will divide the values by these factors to get to the values passed to
% fmincon and multiply by them to get back to the physical values. We do
% this because it seems that fmincon expects all variables to be on roughly
% similar orders of magnitude.
vocr_conv = 1;
phox_conv = 6.25e6;
alpha_conv = 0.01;
conv_vector = [vocr_conv; phox_conv; alpha_conv];

x0 = [1; 6.25e6; 0.04] ./ conv_vector;
lb = [0; 0; 0] ./ conv_vector;
ub = [Inf; Inf; 1] ./ conv_vector;

ss_oh = nan;
ss_tau = nan;

warning('off', 'symbolic:numeric:NumericalInstability')
soln = fmincon(@cost_function, x0, [], [], [], [], lb, ub, [], opts);
warning('on', 'symbolic:numeric:NumericalInstability');

soln = soln .* conv_vector;

ss_vocr = soln(1);
ss_phox = soln(2);
ss_alpha = soln(3);

    function cost = cost_function(x)
        disp([wrf_nox; x])
        if any(x<0)
            error('x cannot be < 0')
        end
        x = x .* conv_vector;
        vocr = x(1);
        phox = x(2);
        alpha = x(3);
        [curr_tau, ~, ~, ~, species] = nox_lifetime(wrf_nox, 'vocr', vocr, 'phox', phox, 'alpha', alpha);
        % compute the cost as a percent error since OH and tau are very
        % different orders of magnitude;
        ss_oh = species.oh;
        ss_tau = curr_tau;
        cost = (abs(reldiff(curr_tau, tau)) + abs(reldiff(species.oh, wrf_oh)))*100;
    end

    function stop = outfun(x,optimvals,state)
        stop = false;
        switch state
            case 'iter'
                x = x .* conv_vector;
                history.x = [history.x; x];
                if show_vals
                    fprintf('VOCR = %.2f, PHOx = %.3g molec. cm^-3, ALPHA = %.2f\n', x(1), x(2), x(3));
                end
        end
    end

end
