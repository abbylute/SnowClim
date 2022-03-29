function[lastpackwater, lastswe, lastpackcc, packsnowdensity,...
    RefrozenWater] = calcEnergyToRefreezing(lastpackwater, lastswe,...
    lastpackcc, lastsnowdepth, WaterDens, LatHeatFreez,...
    RefrozenWater, packsnowdensity)

    % Distribute excess energy to allow refreezing of rain and meltwater 
    % in the pack

    % If packwater > 0 and there is cold content remaining, use it for refreezing
    b = lastpackwater > 0 & lastswe > 0 & lastpackcc < 0 & packsnowdensity < 550;
    % potential energy from refreezing:
    Prf = zeros(size(lastpackwater));
    Prf(b) = WaterDens .* LatHeatFreez .* lastpackwater(b);

    % 1. if available cold content exceeds the potential refreezing
    % energy, then all water is refrozen and cold content is reduced:
    bb = b & -lastpackcc >= Prf;
    lastswe(bb) = lastswe(bb) + lastpackwater(bb);
    RefrozenWater(:,bb) = lastpackwater(bb);
    lastpackwater(bb) = 0;
    lastpackcc(bb) = lastpackcc(bb) + Prf(bb); 
    Prf(bb) = 0;

    % 2. if cold content can't satisfy the refreezing energy, then
    % freeze what you can 
    bb = lastpackwater > 0 & lastswe > 0 & lastpackcc < 0 & -lastpackcc<Prf & packsnowdensity < 550;
    RefrozenWater(:,bb) = -lastpackcc(bb)./(WaterDens .* LatHeatFreez);            
    lastswe(bb) = lastswe(bb) + RefrozenWater(:,bb);
    lastpackwater(bb) = lastpackwater(bb) - RefrozenWater(:,bb);
    lastpackcc(bb) = 0;
    
    % update snow density following refreezing, assuming depth is unchanged.
    packsnowdensity = lastswe .* WaterDens ./ lastsnowdepth;


end