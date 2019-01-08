function nox_lifetime_isopleth(varargin)
%NOX_LIFETIME_ISOPLETH Make an isopleth plot of NOx lifetime
%
%   NOX_LIFETIME_ISOPLETH() Creates an isopleth of NOx lifetime for NOx
%   concentrations between 2e9 and 2e11 molec. cm^-3 and VOCR values
%   between 0.1 and 10.
%
%   Parameters:
%       'phox' - set the P(HOx) value in molec. cm^-3 s^-1. Default is
%       6.25e6.
%
%       'alpha' - set the value of alpha. Default is 0.04.

p = advInputParser;
p.addParameter('phox', 6.25e6);
p.addParameter('alpha', 0.04);
p.parse(varargin{:});
pout = p.Results;

alpha = pout.alpha;
phox = pout.phox;

vocrs = [0.1 1:10];
nox = logspace(-10,-8,10) * 2e19;
taus = nan(length(nox), length(vocrs));

for i_voc = 1:numel(vocrs)
    taus(:,i_voc) = nox_lifetime(nox, 'vocr', vocrs(i_voc), 'alpha', alpha, 'phox', phox);
end

figure; 
nox_mat = repmat(nox(:),1,numel(vocrs));
vocr_mat = repmat(vocrs,numel(nox),1);
contourf(nox_mat, vocr_mat, taus);
cb=colorbar;
cb.Label.String = '\tau (h)';
set(gca,'xscale','log')
xlabel('[NOx]');ylabel('VOC_R')
title(sprintf('\\alpha = %.2f, P(HO_x) = %.2g', alpha, phox));
end

