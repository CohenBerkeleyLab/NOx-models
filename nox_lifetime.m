function [ tau, tau_hno3, tau_ans ] = nox_lifetime( varargin )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

p = inputParser;
p.addOptional('nox', logspace(-11,-7,1000)*2e19, @(x) (isvector(x) && isnumeric(x)));
p.addParameter('phox', 6.25e6, @(x) (isscalar(x) && isnumeric(x)));
p.addParameter('vocr', 5.8,  @(x) (isscalar(x) && isnumeric(x)));
p.addParameter('alpha', 0.04,  @(x) (isscalar(x) && isnumeric(x)));
p.addParameter('no2_no', 4,  @(x) (isscalar(x) && isnumeric(x)));
p.addParameter('k4', 1.1e-11,  @(x) (isscalar(x) && isnumeric(x)));
p.addParameter('k2eff', 8e-12,  @(x) (isscalar(x) && isnumeric(x)));
p.addParameter('k5eff', 5e-12,  @(x) (isscalar(x) && isnumeric(x)));

p.parse(varargin{:});
pout = p.Results;

nox = pout.nox;
phox = pout.phox;
vocr = pout.vocr;
alpha = pout.alpha;
no2_no = pout.no2_no;
k4 = pout.k4;
k2eff = pout.k2eff;
k5eff = pout.k5eff;

args = varargin;
if numel(args)>0 && ~ischar(args{1})
    args = args(2:end);
end

oh = nan(size(nox));
ho2 = nan(size(nox));
ro2 = nan(size(nox));
for a=1:numel(nox)
    [oh(a), ho2(a), ro2(a)] = hox_ss_solver(nox(a), phox, vocr, alpha);
end

T = 298; % kelvin
M = 2e19; % molec. cm^-3

HOME=getenv('HOME');
addpath(fullfile(HOME,'Documents','MATLAB','Rates'));
no = 1/(no2_no + 1) * nox;
no2 = no2_no/(no2_no + 1) * nox;

k_OHNO2 = KOHNO2a(T,M);
k_HO2NO = KNOHO2(T,M);
k_RO2HO2 = 8e-12; % from Paul Romer
k_RO2RO2 = 6.8e-14; % from Paul Romer

tau_hno3 = nox ./ (k_OHNO2 .* oh .* no2 .* 3600).^-1; % convert to hours
% tau_ans = 1 / (alpha * k_RO2+NO * [RO2])
% From Murphy 2006, ACPD, we assume some effective k_RO2_NO (k2eff) that is
% the weighted average of k's for various RO2+NO reactions. [RO2] is
% calculated by assuming:
%  1) RO2 + HO2 is in steady state
%  2) All RO2's go on to produce HO2, thus [RO2] = [HO2]
%
% This means P(RO2) = VOCR*[OH] == L(HO2) = k_HO2+NO * [HO2] * [NO]
% => [HO2] = [RO2] = VOCR * [OH] / (k_HO2+NO * [NO])
%
% This gives us an expression for [RO2], however, we need k_RO2+NO. In
% Murphy, k2eff is the effective rate of NO oxidation by peroxy radicals,
% so
%
%   k2_eff = ( [HO2]*k_HO2+NO + \sum_i (k_RO2_i * [RO2]_i) )/( [HO2 + \sum_i [RO2]_i )
%
% But [HO2] = \sum_i [RO2]_i (let's call this C) so that this becomes:
%
%   k2_eff = ( C*k_HO2+NO + C*k_RO2+NO ) / 2C 
% => 2*k2_eff - k_HO2+NO = k_RO2+NO
k_RO2NO = 2*k2eff - k_HO2NO;
tau_ans = nox ./ (alpha .* k_RO2NO .* ro2 .* no .* 3600);

tau = (1./tau_hno3 + 1./tau_ans).^(-1);

end

