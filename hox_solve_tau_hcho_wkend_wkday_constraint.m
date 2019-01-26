function [wkday_soln, wkend_soln] = hox_solve_tau_hcho_wkend_wkday_constraint(wkday_nox, wkday_tau, wkend_nox, wkend_tau, hcho, alpha, varargin)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
p = advInputParser;
p.addParameter('fmincon_debug', 'none');
p.parse(varargin{:});
pout = p.Results;

E = JLLErrors;

fmincon_debug = pout.fmincon_debug;

vocr_init = 1;
phox_init = 6.25e6;
% In the minimization we'll put P(HOx) in units of ppt/s so that it is a
% similar magnitude to VOCR. It seems that fmincon doesn't do well when the
% values are vastly different orders of magnitude.
phox_convert = 2e19*1e-12;
lb = zeros(3,1);
ub = [100; 1e9; 1e9]; % these values are high enough that no realistic atmosphere should reach them, but they should keep the solver in check

x0 = [vocr_init; phox_init ./ phox_convert; phox_init ./ phox_convert];

last_wkday_tau = nan;
last_wkend_tau = nan;
last_wkday_hcho = nan;
last_wkend_hcho = nan;
history.x = [];

opts = optimoptions('fmincon','Display',fmincon_debug, 'OutputFcn',@outfun);

[fmin_out, ~, fmin_flag] = fmincon(@cost_fxn, x0, [], [], [], [], lb, ub, [], opts);%, -1, 0);
wkday_soln = struct('oh', [], 'ho2', [], 'ro2', [], 'vocr', fmin_out(1),...
    'phox', fmin_out(2) .* phox_convert, 'alpha', alpha, 'tau', [],...
    'last_hcho', last_wkday_hcho, 'last_tau', last_wkday_tau, 'fmincon_flag', fmin_flag);
wkend_soln = struct('oh', [], 'ho2', [], 'ro2', [], 'vocr', fmin_out(1),...
    'phox', fmin_out(3) .* phox_convert, 'alpha', alpha, 'tau', [],...
    'last_hcho', last_wkend_hcho, 'last_tau', last_wkend_tau, 'fmincon_flag', fmin_flag);

[wkday_soln.oh, wkday_soln.ho2, wkday_soln.ro2] = hox_ss_solver(wkday_nox, wkday_soln.phox, wkday_soln.vocr, wkday_soln.alpha);
[wkend_soln.oh, wkend_soln.ho2, wkend_soln.ro2] = hox_ss_solver(wkend_nox, wkend_soln.phox, wkend_soln.vocr, wkend_soln.alpha);
wkday_soln.tau = nox_lifetime(wkday_nox, 'phox', wkday_soln.phox, 'alpha', wkday_soln.alpha, 'vocr', wkday_soln.vocr);
wkend_soln.tau = nox_lifetime(wkend_nox, 'phox', wkend_soln.phox, 'alpha', wkend_soln.alpha, 'vocr', wkend_soln.vocr);

    function stop = outfun(x,optimvals,state)
        stop = false;
        switch state
            case 'iter'
                history.x = [history.x; x];
        end
    end


    function cost = cost_fxn(x)
        vocr = x(1);
        wkday_phox = x(2) .* phox_convert; % back from ppt/s to molec. cm^-3 s^-1
        wkend_phox = x(3) .* phox_convert;
        
        last_wkday_tau = nox_lifetime(wkday_nox, 'phox', wkday_phox, 'alpha', alpha, 'vocr', vocr);
        last_wkend_tau = nox_lifetime(wkend_nox, 'phox', wkend_phox, 'alpha', alpha, 'vocr', vocr);
        lifetime_cost = (last_wkday_tau - wkday_tau).^2 + (last_wkend_tau - wkend_tau).^2;
        
        % the second part of the cost is the difference between the sat +
        % profile inferred HCHO and the steady state model. we convert from
        % molec. cm^-3 to ~ppb so that it is a similar magnitude to the
        % lifetime cost.
        last_wkday_hcho = hcho_steady_state(vocr, wkday_phox);
        last_wkend_hcho = hcho_steady_state(vocr, wkend_phox);
        hcho_cost =  (last_wkday_hcho - hcho).^2 + (last_wkend_hcho - hcho).^2;
        cost = lifetime_cost + hcho_cost .* 1e9.^2 ./ 2e19.^2;
    end

end

