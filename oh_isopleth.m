function oh_isopleth(varargin)
%OH_ISOPLETH Make an isopleth plot of OH concentration.
%
%   OH_ISOPLETH() Creates an isopleth of OH concentration for NOx
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
oh = nan(length(nox), length(vocrs));

for i_voc = 1:numel(vocrs)
    for i_nox = 1:numel(nox)
        oh(i_nox,i_voc) = hox_ss_solver(nox(i_nox), phox, vocrs(i_voc), alpha);
    end
end

figure; 
nox_mat = repmat(nox(:),1,numel(vocrs));
vocr_mat = repmat(vocrs,numel(nox),1);
contourf(nox_mat, vocr_mat, log10(oh));
cb=colorbar;
cbticks = floor(cb.Ticks(1)):ceil(cb.Ticks(end));
cb.Ticks = cbticks;
cb.TickLabels = sprintfmulti('10^{%d}', cbticks);
caxis([cbticks(1), cbticks(end)]);

cb.Label.String = '[OH] (molec. cm^{-3})';
set(gca,'xscale','log')
xlabel('[NOx]');ylabel('VOC_R')
colormap jet

end

