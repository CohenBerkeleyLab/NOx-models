function [ C ] = advection_model(  )
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


dt = 1; % second
dx = 1e5; % cm
dy = 1e5;
dz = 5e4;

P = 101325; % Pa
T = 298; % K

u = 300; % cm/s 

nsteps = 86400;
nboxes = 100;

doemiss = true;
emiss = 250 / 1e10 / 3600 * 6.02e23; % 250 mol/km^2*hr --> molec/cm^2*s
emiss = emiss * 0.01;

doloss = true;
k_loss = 1e-11;

C = nan(nboxes, nsteps+1);

Nair = P ./ (T .* 8.314) .* 6.022e23 ./ 1e6;

C(1, 1) = 1e-9 * Nair;
C(2:end, 1) = 0;

for t=1:nsteps
    dC_advect = advect(C(:,t)) .* dt;
    if doemiss
        dC_emiss = emiss .* dt ./ dz;
    else
        dC_emiss = 0;
    end
    if doloss
        dC_chem = -k_loss .* C(1,t) .* dt;
    else 
        dC_chem = 0;
    end
    
    C(:,t+1) = C(:,t) + dC_advect + dC_emiss + dC_chem;
end

    function dCdt = advect(Cvec)
        % dC/dt = -u * dC/dx
        dCdt = nan(size(Cvec));
        dCdt(1) = -u .* Cvec(1) ./ dx;
        for x=2:nboxes
            dCdt(x) = -u .* (Cvec(x) - Cvec(x-1)) ./ dx;
        end
    end

end

