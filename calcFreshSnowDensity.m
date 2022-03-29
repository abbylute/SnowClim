function newsnowdensity = calcFreshSnowDensity(airtemp)

% calculate the density of fresh snowfall as a function of air temperature
% based on eqn 16 from Essery et al., (2013), which is from Anderson
% (1976).
% Constants are adopted from Table 4 in Essery et al., (2013), which are
% taken from Oleson et al., (2004).

% airtemp should be in degrees C

% define constants:
df = 1.7;   % 1/K
ef = 15;    % K
pmin = 50;  % kg/m3 minimum snow density

% temperatures must be > ef, otherwise we get imaginary numbers
airtemp(airtemp < -ef) = ef;

newsnowdensity = pmin + max(df * (airtemp + ef).^(3/2), 0);
