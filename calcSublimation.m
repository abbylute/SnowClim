function [sublimation, condensation, SnowWaterEq, SnowDepth, lastpackcc,...
    SnowDensity, lastpackwater] = calcSublimation(E, SnowWaterEq, ...
    SnowDepth, SnowDensity, lastsnowtemp, lastpackcc, SnowDensDefault,...
    sublimation, condensation, lastpackwater, WaterDens)

    % latent heat of sublimation
    lambdaS = 2834.1 - 0.29*lastsnowtemp - 0.004*lastsnowtemp.^2; 

    % calculate sublimation and evaporation
    Sublimation = zeros(1,length(E));
    Evaporation = zeros(1,length(E),1);
    Sublimation(lastsnowtemp<0) = -E(lastsnowtemp<0)./WaterDens; % E has units of kg/m2/s
    Evaporation(lastsnowtemp==0) = -E(lastsnowtemp==0)./WaterDens; % E has units of kg/m2/s

    % Update SWE, depth, density, cold content
        f=find(SnowWaterEq>Sublimation); % includes possibility of condensation (ie negative sublimation)
        f1=find(SnowWaterEq<=Sublimation);

        % no limit
        initialSWE = SnowWaterEq(f);
        SnowWaterEq(f) = SnowWaterEq(f)-Sublimation(f);
        SnowDepth(f) = SnowWaterEq(f)./SnowDensity(f)*WaterDens;
        %f = find(SnowWaterEq>Sublimation & Sublimation>0); % only update cc for sublimation, not condensation
        %lastpackcc(f) = lastpackcc(f).* (SnowWaterEq(f)./initialSWE);

        % Sublimation eats all of the snow
        Sublimation(f1) = SnowWaterEq(f1);
        SnowWaterEq(f1) = 0;
        SnowDepth(f1) = 0;
        lastpackcc(f1) = 0;
        SnowDensity(f1) = SnowDensDefault;    
       
    % Output sublimation and condensation    
    b = Sublimation > 0;
    sublimation(b) = sublimation(b) + Sublimation(b);
    condensation(~b) = condensation(~b) + Sublimation(~b);    
    
    % update packwater
    lastpackwater = max(0, lastpackwater - Evaporation);
end