function [ oh_nd, ho2_nd, ro2_nd ] = hox_ss_solver( nox, phox, vocr, alpha )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

HOME=getenv('HOME');
addpath(fullfile(HOME,'Documents','MATLAB','Rates'));

T = 298;
M = 2e19;

k_RO2NO = 8e-12; % approximate guess
k_HO2NO = KNOHO2(T,M);
k_HO2HO2 = KHO2self(T,M, 0.01*M); % assume 1% humidity
k_RO2HO2 = 8e-12;
k_RO2RO2 = 6.8e-14;
k_OHNO2 = KOHNO2a(T,M);

no = 0.2*nox;
no2 = 0.8*nox;

syms ro2 ho2 oh

eqns = [ho2 == (k_RO2NO * ro2 * no) / (k_HO2NO * no + k_HO2HO2 * ho2 + k_RO2HO2 * ro2),...
        ro2 == (vocr * oh) / (k_RO2NO * no + k_RO2HO2 * ho2 + k_RO2RO2 * ro2),...
        phox == k_OHNO2 * oh * no2 + alpha * k_RO2NO * ro2 * no + k_RO2HO2 * ro2 * ho2 * k_RO2RO2 * ro2.^2 + k_HO2HO2 * ho2.^2];

    
S = vpasolve(eqns, [ho2, ro2, oh]);
ho2_nd = double(S.ho2);
ro2_nd = double(S.ro2);
oh_nd = double(S.oh);
test_mat = [ho2_nd, ro2_nd, oh_nd];
rr = false(size(test_mat));
pp = real(test_mat) > 0;
for a=1:numel(test_mat)
    rr(a) = isreal(test_mat(a));
end

xx = all(pp & rr,2);
ho2_nd = ho2_nd(xx);
ro2_nd = ro2_nd(xx);
oh_nd = oh_nd(xx);

end

