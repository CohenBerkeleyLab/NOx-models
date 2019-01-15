function [oh, ho2, ro2, soln] = hox_solve_tau_hcho_constraint(nox, tau, hcho, alpha, varargin)
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
lb = zeros(2,1);
ub = Inf(2,1);

x0 = [vocr_init; phox_init ./ phox_convert];

last_tau = nan;
last_hcho = nan;
history.x = [];
opts = optimoptions('fmincon','Display',fmincon_debug,'OutputFcn',@outfun);

fmin_out = fmincon(@cost_fxn, x0, [], [], [], [], lb, ub, [], opts);%, -1, 0);
soln = struct('vocr', fmin_out(1), 'phox', fmin_out(2) .* phox_convert,...
    'alpha', alpha, 'tau', [], 'last_hcho', last_hcho, 'last_tau', last_tau);

[oh, ho2, ro2] = hox_ss_solver(nox, soln.phox, soln.vocr, soln.alpha);
soln.tau = nox_lifetime(nox, 'phox', soln.phox, 'alpha', soln.alpha, 'vocr', soln.vocr);

    function stop = outfun(x,optimvals,state)
        stop = false;
        switch state
            case 'iter'
                history.x = [history.x; x];
        end
    end


    function cost = cost_fxn(x)
        vocr = x(1);
        phox = x(2) .* phox_convert; % back from ppt/s to molec. cm^-3 s^-1
        last_tau = nox_lifetime(nox, 'phox', phox, 'alpha', alpha, 'vocr', vocr);
        lifetime_cost = last_tau - tau;
        % the second part of the cost is the difference between the sat +
        % profile inferred HCHO and the steady state model. we convert from
        % molec. cm^-3 to ~ppt so that it is a similar magnitude to the
        % lifetime cost.
        last_hcho = hcho_steady_state(vocr, phox);
        hcho_cost =  (last_hcho - hcho) ./ 2e19 .* 1e12;
        cost = lifetime_cost .^ 2 + hcho_cost .^ 2;
    end

    function hcho = hcho_steady_state(vocr, phox)
        % This steady state model is from Valin et al. JGR 2016 (doi:
        % 10.1002/2015JD024012). The values for alpha effective, j_hcho,
        % and k_hcho are from that paper as well. Note that alpha effective
        % is not the RO2 branching ratio but the fractional production of
        % HCHO from all reactions that consume OH.
        
        % phox = molec. cm-3 s^-1
        % vocr = s^-1
        hcho_alpha_eff = 0.3; % unitless
        j_hcho = 8.3e-5; % s^-1
        k_hcho = 8.3e-12; % cm^3 molec.^-1 s^-1
        
        hcho = (hcho_alpha_eff .* phox) ./ (j_hcho + k_hcho .* phox ./ vocr);
    end

end

