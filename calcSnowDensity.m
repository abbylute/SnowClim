function packsnowdensity = calcSnowDensity(lastswe, lastsnowtemp,...
    packsnowdensity, WaterDens, sec_in_ts)

% calculate new snow density following compaction as a function of swe and
% snow temperature,
% based on equations 11 and 12 in Essery et al., (2013), which are based on
% Anderson (1976). Constants are adopted from Table 3 in Essery et al.,
% (2013), which are based on Boone (2002) and the ISBA-ES snow model.

% lastsnowtemp should be in degrees C

% define constants:
c1 = 2.8E-6; % 1/s
c2 = 0.042;   % 1/K
c3 = 0.046;   % m3/kg
c4 = 0.081;   % 1/K
c5 = 0.018;   % m3/kg
p0 = 150;     % kg/m3
n0 = 3.7E7;   % kg/m/s can vary greatly in the literature from 3.6E6 to 3.7E7
g  = 9.81;    % m/s2


% use average snow mass between ground surface and snow surface
% (half of total mass)
mass = (lastswe/2) * WaterDens;

% calculate viscosity, nu:
nu = n0 * exp(c4*-lastsnowtemp + c5*packsnowdensity);

% calculate change in snow density:
delta_density = packsnowdensity .* ((mass*g)./nu + c1*exp(-c2*-lastsnowtemp - c3*max(0, packsnowdensity-p0)));

% convert delta_density from kg m-3 s-1 to appropriate timestep
delta_density = delta_density*sec_in_ts;

% calculate new snow density:
packsnowdensity = packsnowdensity + delta_density;








