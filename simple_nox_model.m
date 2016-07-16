function [ c ] = simple_nox_model( )
%SIMPLE_NOX_MODEL Very simple 1 box model for comparison with that from WRF
%   Basically I need a very simple NOx model to compare to my WRF parsing
%   code. Right now, just want to simulate:
%
%       NO + O3 --> NO2
%       NO2 + hv --> O3
%
%   will add more reactions later

nt = 10800;
dt = 1; % seconds

no = nan(1,nt);
no2 = nan(1,nt);
o3 = nan(1,nt);

nair = 2e19; % molec./cm^3
T = 298; % K

no(1) = 1e-9*nair;
no2(1) = 1e-9*nair;
o3(1) = 40e-9*nair;

HOMEDIR=getenv('HOME');
addpath(fullfile(HOMEDIR,'Documents','MATLAB','Rates'));

% Get photolysis rates for noon
phot_date = '2013-06-01';
phot_lon = -84;
phot_lat = 34;

jNO2 = call_tuv('Pj_no2',phot_date,12,phot_lon,phot_lat,false);

for t=1:nt-1
    dno = -KNOO3(T,nair)*no(t)*o3(t) + no2(t)*jNO2;
    dno2 = KNOO3(T,nair)*no(t)*o3(t) + -no2(t)*jNO2;
    do3 = -KNOO3(T,nair)*no(t)*o3(t) + no2(t)*jNO2;
    
    no(t+1) = no(t)+dno*dt;
    no2(t+1) = no2(t)+dno2*dt;
    o3(t+1) = o3(t)+do3*dt;
end

c(1,:) = no./nair;
c(2,:) = no2./nair;
c(3,:) = o3./nair;

end

