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