function packsnowdensity = calcSnowDensityAfterSnow(lastswe, newswe,...
    packsnowdensity, newsnowdensity)

% Calculate the new density of the snowpack following fresh snowfall.
% Based on Essery et al., (2013) eqn 17.


packsnowdensity = (lastswe + newswe)./...
    ((lastswe)./packsnowdensity + (newswe)./newsnowdensity);
