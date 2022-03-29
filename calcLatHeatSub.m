function[LatHeatSub] = calcLatHeatSub(lastsnowtemp)

% latent heat of sublimation (kJ/kg)
        LatHeatSub = 2834.1 - 0.29*lastsnowtemp - 0.004*lastsnowtemp.^2;
end
