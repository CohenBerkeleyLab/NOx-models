function oh_vs_lifetime(nox, vocr, varargin)
%OH_VS_LIFETIME Plot how [OH] and NOx lifetime vary against each other
%   OH_VS_LIFETIME( NOx, VOCR ) Given a NOx concentration in molec. cm^-3
%   and a VOCR value in s^-1, this will vary the NOx concentration by plus
%   or minus one order of magnitude and the VOCR by +/- 5 (staying above 0)
%   and plot the steady state OH concentration and NOx lifetime for each
%   point.

p10 = log10(nox);
nox_vec = unique([logspace(p10-1, p10+1, 10), nox]);
vocr_vec = (vocr-5):(vocr+5);
vocr_vec = vocr_vec(vocr_vec>0);

[tau_nox, ~, ~, ~, species_nox] = nox_lifetime(nox_vec, 'vocr', vocr, varargin{:});
oh_nox = species_nox.oh;

tau_vocr = nan(size(vocr_vec));
oh_vocr = nan(size(vocr_vec));
for i=1:numel(vocr_vec)
    [tau_vocr(i),~,~,~,spec] = nox_lifetime(nox, 'vocr', vocr_vec(i), varargin{:});
    oh_vocr(i) = spec.oh;
end

figure;
subplot_stretch(1,2);
subplot(1,2,1);
[~,i_nox] = min(abs(nox_vec - nox));
scatter(tau_nox, oh_nox, [], nox_vec,'filled');
cb=colorbar;
cb.Label.String = '[NO_x]';
line(tau_nox(i_nox), oh_nox(i_nox), 'marker','o', 'color', 'k', 'markersize',10,'linewidth',2,'linestyle','none');
xlabel('\tau (h)');
ylabel('[OH] (molec. cm^{-3})')
title('Varied NOx')


subplot(1,2,2);
[~,i_vocr] = min(abs(vocr_vec - vocr));
scatter(tau_vocr, oh_vocr, [], vocr_vec,'filled');
cb=colorbar;
cb.Label.String = 'VOC_R';
line(tau_vocr(i_vocr), oh_vocr(i_vocr), 'marker','o', 'color', 'k', 'markersize',10,'linewidth',2,'linestyle','none');
xlabel('\tau (h)');
ylabel('[OH] (molec. cm^{-3})')
title('Varied VOCR')
end

