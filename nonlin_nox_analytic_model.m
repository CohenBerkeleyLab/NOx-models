function [ oh ] = nonlin_nox_analytic_model( nox, varargin )
%OH = NONLIN_NOX_ANALYTIC_MODEL( NOX, ... ) Analytical model for nonlinear OH-NOx chemistry
%   Source: Murphy et al. ACPD, 2006, 6, 11971
%   This function will return the OH concentrations corresponding to the
%   solution of the analytical model given in Murphy et al 2006 for the
%   vector of NOx concentrations input.  By default it will use the values
%   given in Table 4 of Murphy for the various parameters, but each one can
%   be overridden by the following parameter arguments:
%
%       'phox' = production of HOx, default = 1.25 ppb/hr = 6.25e6 molec.
%           cm^-3 s^-1 assuming number density of air is 2e19 molec. cm^-3
%       'vocr' = k1*[VOC] i.e. VOC reactivity with OH, default = 5.8 s^-1
%       'alpha' = RONO2 branching ratio, default = 0.04.
%       'no2/no' = NO2 to NO ratio, default = 4.
%       'k4' = rate constant for reaction of OH + NO2 --> HNO3, default =
%           1.1e-11 s^-1 cm^3 molec.^-1
%       'k2eff' = effective reaction rate of NO with RO2, default = 8e-12
%           s^-1 cm^3 molec.^-1
%       'k5eff' = effective rate of peroxy (RO2, HO2) self reaction,
%           default = 5e-12 s^-1 cm^3 molec.^-1
%
%   Obviously given the unit of the rate constants, concentrations of NOx
%   should be input in molec./cm^3.

p = inputParser;
p.addParameter('phox', 6.25e6, @(x) (isscalar(x) && isnumeric(x)));
p.addParameter('vocr', 5.8,  @(x) (isscalar(x) && isnumeric(x)));
p.addParameter('alpha', 0.04,  @(x) (isscalar(x) && isnumeric(x)));
p.addParameter('no2_no', 4,  @(x) (isscalar(x) && isnumeric(x)));
p.addParameter('k4', 1.1e-11,  @(x) (isscalar(x) && isnumeric(x)));
p.addParameter('k2eff', 8e-12,  @(x) (isscalar(x) && isnumeric(x)));
p.addParameter('k5eff', 5e-12,  @(x) (isscalar(x) && isnumeric(x)));

p.parse(varargin{:});
pout = p.Results;

phox = pout.phox;
vocr = pout.vocr;
alpha = pout.alpha;
no2_no = pout.no2_no;
k4 = pout.k4;
k2eff = pout.k2eff;
k5eff = pout.k5eff;

% After all that, it's really just one long equation. Well, after we go
% ahead and compute the NO2 and NO concentrations, but that's to simple to
% mention.

no = 1/(no2_no + 1) * nox;
no2 = no2_no/(no2_no + 1) * nox;

% Okay, it's actually just the quadratic formula, so to hopefully eliminate
% typos, I'll just define a, b, and c and use those.

a = 6 .* k5eff .* (vocr ./ (k2eff .* no)).^2;
b = k4 .* no2 + alpha .* vocr;
c = -phox;

oh = ( -b + sqrt( b.^2 - 4 .* a .* c )) ./ (2 .* a);


end

