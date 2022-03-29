function[SnowMelt, MeltEnergy, lastpackwater, lastswe, lastsnowdepth]...
    = calcEnergyToMelt(lastswe, lastsnowdepth, packsnowdensity, lastenergy,...
    lastpackwater, LatHeatFreez, WaterDens, SnowMelt, MeltEnergy)

    % Distribute excess energy to snowmelt

    b = lastswe>0 & lastenergy>0;       
    potmelt = lastenergy(b) ./ (LatHeatFreez * WaterDens); 
    melt = zeros(1,length(lastpackwater));
    melt(b) = min([lastswe(b); potmelt],[],1);
    melt(melt<0) = 0;
    SnowMelt(:,b) = melt(b); % add to the monthly melt
    lastsnowmelt=melt;
    meltenergy = melt(b) .* LatHeatFreez .* WaterDens;
    lastenergy(b) = lastenergy(b) - meltenergy;
    MeltEnergy(:,b) = meltenergy;

    % Update water in snowpack
    lastpackwater = lastpackwater + melt;
    
    
    % Update snow
    lastswe = lastswe - melt;
    b = lastswe > 0;
    lastsnowdepth(b) = lastswe(b) * WaterDens ./ packsnowdensity(b);
    lastsnowdepth(~b) = 0;

end