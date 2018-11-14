function nox_lifetime_isopleth()
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
vocrs = [0.1 1:10];
nox = logspace(-10,-8,10) * 2e19;
taus = nan(length(nox), length(vocrs));

for i_voc = 1:numel(vocrs)
    taus(:,i_voc) = nox_lifetime(nox, 'vocr', vocrs(i_voc));
end

figure; 
nox_mat = repmat(nox(:),1,numel(vocrs));
vocr_mat = repmat(vocrs,numel(nox),1);
contourf(nox_mat, vocr_mat, taus);
colorbar;
set(gca,'xscale','log')
xlabel('[NOx]');ylabel('VOC_R')

end

