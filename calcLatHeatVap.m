function[LatHeatVap] = calcLatHeatVap(lastsnowtemp)       

% latent heat of vaporization (kJ/kg)
        LatHeatVap = 2500.8 - 2.36*lastsnowtemp + 0.016*lastsnowtemp.^2 + ...
            0.00006*lastsnowtemp.^3;
end
