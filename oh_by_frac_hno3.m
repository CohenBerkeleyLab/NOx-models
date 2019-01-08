function fig = oh_by_frac_hno3(varargin)
%OH_BY_FRAC_HNO3 Plot the relationship between [NOx], [OH], and frac. loss to HNO3
%
%   OH_BY_FRAC_HNO3() plot with default values for NOx, VOCR, PHOx, and
%   alpha.
%
%   Parameters:
%
%       'nox' - vector of NOx concentrations to plot for in molec. cm^-3.
%
%       'vocr' - vector of VOCR values to plot (each one is a series) in
%       s^-1.
%
%       'phox' - scalar value of P(HOx) to use for the steady state
%       calculations in molec. cm^-3 s^-1.
%
%       'alpha' - scalar value of alpha to use for the steady state
%       calculations.
%
%       'fracpos' - controls where fraction of NOx lost to HNO3 is shown.
%       Default is 'color', i.e. is plotted as the color. May also be 'x',
%       to plot it on the x axis and [NOx] on the colorbar.

p = advInputParser;
p.addOptional('nox', logspace(-1,2,30)*2e10);
p.addOptional('vocr', [1, 5, 10]);
p.addParameter('phox', 6.25e6);
p.addParameter('alpha', 0.04);
p.addParameter('fracpos', 'color');

p.parse(varargin{:});
pout = p.Results;

E = JLLErrors;

nox = pout.nox;
vocr = pout.vocr;
phox = pout.phox;
alpha = pout.alpha;
frac_pos = pout.fracpos;

n_nox = numel(nox);
n_vocr = numel(vocr);

markers = {'o','^','v', 'p'};

l = gobjects(n_vocr,1);
fig = figure;
hold on

for i_vocr = 1:n_vocr
    [~, tau_hno3, tau_ans, ~, species] = nox_lifetime(nox, 'vocr', vocr(i_vocr), 'phox', phox, 'alpha', alpha);
    k_hno3 = 1 ./ tau_hno3;
    k_ans = 1 ./ tau_ans;
    
    oh_conc = species.oh;
    frac_hno3 = k_hno3 ./ (k_hno3 + k_ans);
    
    switch lower(frac_pos)
        case 'color'
            l(i_vocr) = scatter(nox, oh_conc, [], frac_hno3, 'filled', 'marker', markers{i_vocr});
            xstr = '[NO_x] (molec. cm^{-3})';
            set(gca, 'xscale', 'log')
            cblabel = 'Frac. loss to HNO_3';
        case 'x'
            l(i_vocr) = scatter(frac_hno3, oh_conc, [], log10(nox), 'filled', 'marker', markers{i_vocr});
            xstr = 'Frac. loss to HNO_3';
            cblabel = 'log10([NO_x]) (molec. cm^{-3})';
        otherwise
            E.notimplemented('colorvar = %s', frac_pos);
    end
end

xlabel(xstr)
ylabel('[OH] (molec. cm^{-3})')
cb = colorbar;
cb.Label.String = cblabel;
legstr = sprintfmulti('VOC_R = %.1f', vocr);
legend(l, legstr, 'location', 'best');

end

