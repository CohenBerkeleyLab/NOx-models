function oh_isopleth()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
vocrs = [0.1 1:10];
nox = logspace(-10,-8,10) * 2e19;
oh = nan(length(nox), length(vocrs));
phox = 6.25e6;
alpha = 0.04;

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

