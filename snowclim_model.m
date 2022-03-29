function [SnowMelt, SnowWaterEq, SFE, SnowDepth, SnowDensity,...
    Sublimation, Condensation, SnowTemp, MeltEnergy, Energy, Albedo,...
    SnowYN, RaininSnow, Runoff, RefrozenWater, PackWater, LW_down, LW_up,...
    SW_down, SW_up ,Q_latent, Q_sensible, Q_precip, PackCC, CCenergy,...
    CCsnowfall]...
    = snowclim_model(lat, lrad, tavg, ppt, solar, tdmean,...
    vs, relhum, psfc, huss, parameterfilename)


disp('running SnowClim model...')

%--- Input Data ---%
    % lat                   # latitudes for points to model at (1 x space)
    % lrad                  # downward longwave radiation (kJ/m2/hr) (time x space)
    % tavg                  # average air temperature (C) (time x space)
    % ppt                   # precipitation (m) (time x space)
    % solar                 # downward shortwave radiation (kJ/m2/hr) (time x space)
    % tdmean                # dewpoint temperature (C) (time x space)
    % vs                    # windspeed (m/s) (time x space)
    % relhum                # relative humidity (%) (time x space)
    % psfc                  # air pressure (hPa or mb) (time x space)
    % huss                  # specific humidity (kg/kg) (time x space)
    
    
%--- Parameters ---%
    load(parameterfilename, 'S') 
    cal = S.cal;

    
%--- Constants ---%
    WaterDens = 1000;       % density of water (kg/m3)
    LatHeatFreez = 333.3;   % latent heat of freezing (kJ/kg) 
    Ci = 2.117;             % heat Capacity of Snow (kJ/kg/C)
	Cw = 4.2;               % heat Capacity of Water (kJ/kg/C)
    Ca = 1.005;             % heat Capacity of Air (kJ/kg/C)
    k = 0.41;               % von Karman's constant
    
    
%--- Allocate space ---%
    % dimensions of output (time, space)
    outdim = size(ppt);
    
    % mass outputs:
    SnowDepth     = single(NaN(outdim)); %  snow depth (m)
    SnowWaterEq   = single(NaN(outdim)); %  snow water equivalent (m)
    SnowMelt 	  = single(NaN(outdim)); %  snow melt (m)
    Sublimation   = single(NaN(outdim)); %  snow sublimation (m)
    Condensation  = single(NaN(outdim)); %  snow condensation (m)
    Runoff        = single(NaN(outdim)); %  runoff from snowpack (m)
    RaininSnow    = single(NaN(outdim)); %  rain added to snowpack (m)
    RefrozenWater = single(NaN(outdim)); %  liquid water that is refrozen in the snowpack (m)
    PackWater     = single(NaN(outdim)); %  liquid water in the snowpack (m)
    SnowYN        = single(NaN(outdim)); %  snow cover binary
    
    % values below are the average value on days with snow on the ground
    SnowTemp 	  = single(NaN(outdim)); % temperature of snow (C)
    Albedo 		  = single(NaN(outdim)); % snow surface albedo
    SnowDensity   = single(NaN(outdim)); % snowpack density (kg/m3)

    % energy outputs:
    Q_sensible    = single(NaN(outdim)); % sensible heat flux (kJ/m2/timestep)
    Q_latent      = single(NaN(outdim)); % latent heat flux (kJ/m2/timestep)
    Q_precip      = single(NaN(outdim)); % precipitation heat flux (kJ/m2/timestep)
    LW_up         = single(NaN(outdim)); % upward longwave radiation from the snow surface (kJ/m2/timestep)
    LW_down       = single(NaN(outdim)); % downward longwave radiation to the snow surface (kJ/m2/timestep)
    SW_up         = single(NaN(outdim)); % upward shortwave radiation from the snow surface (kJ/m2/timestep)
    SW_down       = single(NaN(outdim)); % downward shortwave radiation to the snow surface (kJ/m2/timestep)
	Energy 	      = single(NaN(outdim)); % net energy to the snowpack (kJ/m2/timestep)
	MeltEnergy 	  = single(NaN(outdim)); % energy used for melting snow (kJ/m2/timestep)
    PackCC        = single(NaN(outdim)); % snowpack cold content (kJ/m2/timestep)
    CCsnowfall    = single(NaN(outdim)); % cold content contributed by snowfall (kJ/m2/timestep)
    CCenergy      = single(NaN(outdim)); % cold content contributed by changes in net energy (kJ/m2/timestep)
       
        
%---Converted Inputs ---%
    % number of seconds in each time step
    sec_in_ts = S.hours_in_ts*3600;

    % calculate % of ppt as snow/rain
    passnow = calcPhase(tavg+273.15, ones(size(tavg)), relhum, false);
    % depth of rain (m)
    R_m=ppt.*(1-passnow); 

    % snowfall water equivalent (m)
    SFE = ppt.*passnow; 
    
    % new snowfall less than 0.1mm/hr is set to 0
    a = SFE < (.0001*S.hours_in_ts); 
    SFE(a) = 0; 
    R_m(a) = ppt(a);
    
    % density of new snowfall
    newsnowdensity = calcFreshSnowDensity(tavg); % (kg/m3)
    
    
%--- For each time step ---%
for i = 1:size(cal,1)
    
    %--- reset to 0 snow at specified time of year ---%
    if (i==1 || (cal(i,2)==S.snowoff_month && cal(i,3)==S.snowoff_day)) 
        lastalbedo = single(ones(1,length(lat))*S.ground_albedo);
        lastpacktemp = single(zeros(1,length(lat)));
        lastswe = single(zeros(1,length(lat)));
        lastsnowdepth = single(zeros(1,length(lat)));
        packsnowdensity = single(ones(1,length(lat))*S.snow_dens_default);
        snowage = single(zeros(1,length(lat)));
        lastpackcc = single(zeros(1,length(lat)));
        lastpackwater = single(zeros(1,length(lat)));
    end
    
    
    %--- new mass inputs ---%
    newswe = SFE(i,:);
    newrain = R_m(i,:);
    newsnowdepth = SFE(i,:) .* WaterDens ./ newsnowdensity(i,:);
    newsnowdens = newsnowdensity(i,:);

    
    %--- Calculate snow temperature and cold contents ---%
    newsnowtemp = single(zeros(1,length(lat)));
    f = newswe > 0;
    newsnowtemp(f) = min(0, tdmean(i,f)); 
    
    snowfallcc = WaterDens .* Ci .* newswe .* newsnowtemp;
    lastpackcc = lastpackcc + snowfallcc;
    CCsnowfall(i,:) = snowfallcc;
    lastpacktemp(f) = lastpackcc(f) ./ (WaterDens .* Ci .* (lastswe(f)+newswe(f)));        


    % If there is snow on the ground, run model 
    a = newswe + lastswe > 0;
    if sum(a)>0             
        SnowYN(i,:) = a;
       
 
        %--- Set snow surface temperature ---%
        lastsnowtemp = min(tdmean(i,:)+S.Ts_add,0);
        SnowTemp(i,a) = lastsnowtemp(a);


        %--- Update snowpack after new snowfall ---%
        packsnowdensity(a) = single(calcSnowDensityAfterSnow(lastswe(a),...
            newswe(a), packsnowdensity(a), newsnowdens(a)));
        lastswe = lastswe + newswe;
        lastsnowdepth = lastsnowdepth + newsnowdepth;
        
        
        %--- Calculate pack density after compaction ---%
        packsnowdensity(a) = calcSnowDensity(lastswe(a), lastpacktemp(a),...
            packsnowdensity(a), WaterDens, sec_in_ts);
        lastsnowdepth(a) = lastswe(a).*WaterDens./packsnowdensity(a); 


        %--- Update snowpack liquid water content ---%
        previouspackwater = lastpackwater(a);
        lastpackwater(a) = lastpackwater(a) + newrain(a);
        [Runoff(i,:), lastpackwater] = updatePackWater(a, lastpackwater,...
            lastsnowdepth, S.lw_max, Runoff(i,:), sec_in_ts);
        RaininSnow(i,a) = max(lastpackwater(a) - previouspackwater, 0);


        %--- Calculate albedo ---%
        [lastalbedo, snowage] = calcAlbedo(S.albedo_option,...
            S.ground_albedo, S.max_albedo, lastalbedo,newsnowdepth,...
            lastsnowdepth, newswe, lastswe, lastsnowtemp, WaterDens, lat,...
            cal(i,2), cal(i,3), snowage, lastpackcc, sec_in_ts);
        Albedo(i,a) = lastalbedo(a);


        %--- Calculate turbulent heat fluxes (kJ/m2/timestep) ---%
        H = single(zeros(size(lat)));
        E = single(zeros(size(lat)));
        EV = single(zeros(size(lat)));
        [H(a), E(a), EV(a)] = calcTurbulentFluxes(S.stability, S.windHt,...
            S.z_0, S.tempHt, S.z_h, k, vs(i,a), lastsnowtemp(a),...
            tavg(i,a), Ca, psfc(i,a), huss(i,a), S.E0_value, S.E0_app,...
            S.E0_stable, sec_in_ts);


        %--- Rain heat flux into snowpack (kJ/m2/timestep) ---%
        P = single(zeros(size(lastsnowtemp)));
        P(a) = Cw .* WaterDens .* max(0,tdmean(i,a)) .* newrain(a);
        Q_precip(i,:) = P;


        %--- Net downward solar flux at surface (kJ/m2/timestep) ---%
        Sup = single(zeros(size(lat)));
        Sup(a) = solar(i,a) .* lastalbedo(a);
        Sdn = single(zeros(size(lat)));
        Sdn(a) = solar(i,a);


        %--- Longwave flux up from snow surface (kJ/m2/timestep) ---%
        Ldn = single(zeros(size(lat)));
        Ldn(a) = lrad(i,a);
        Lt = single(zeros(size(lat)));
        Lt(a) = calcLongwave(S.snow_emis, lastsnowtemp(a), lrad(i,a),...
            sec_in_ts);


        %--- Ground heat flux (kJ/m2/timestep)---%
        Gf = single(zeros(size(lat)));
        Gf(a) = S.G * sec_in_ts;


        %--- Downward net energy flux into snow surface (kJ/m2/timestep) ---%
        lastenergy = Sdn - Sup + Ldn - Lt + H + EV + Gf + P;
        Energy(i,:) = lastenergy;

        
        %--- Apply cold content tax ---%
        lastenergy = nanmean(Energy(max(1,(i-S.smooth_hr+1)):i,:),1);
        tax = (lastpackcc - S.Tstart)./(S.Tadd).*S.maxtax;
        % limit tax to be >= 0 and <= maxtax
        tax = max(0, tax); 
        tax = min(S.maxtax, tax);
        % apply tax to negative energy only
        n = lastenergy<0;
        lastenergy(n) = lastenergy(n) .* (1-tax(n));
        
        
        %--- Distribute energy ---%

            %--- 1. Energy goes to cold content first ---%
            [lastpackcc, lastenergy, CCenergy(i,:)]...
                = calcEnergyToCC(lastpackcc, lastenergy, CCenergy(i,:));
            b = find(lastswe > 0);            
            lastpacktemp(b) = lastpackcc(b)./ (WaterDens .* Ci .* lastswe(b));

            
            %--- Apply temperature instability correction ---%
            thres = sec_in_ts/3600 * 0.015; % 15mm for each hour in the time step
            f = find(lastswe < thres & lastswe > 0 & lastpacktemp < (tavg(i,:)));
            lastpacktemp(f) = min(0, tavg(i,f));
            lastpackcc(f) = WaterDens .* Ci .* lastswe(f) .* lastpacktemp(f);

            
            %--- 2. Energy goes to refreezing second ---%
            [lastpackwater, lastswe, lastpackcc, packsnowdensity,...
                RefrozenWater(i,:)] = calcEnergyToRefreezing(lastpackwater,...
                lastswe, lastpackcc, lastsnowdepth, WaterDens, LatHeatFreez,...
                RefrozenWater(i,:), packsnowdensity);
            

            %--- 3. Energy goes to melt third ---% 
            [SnowMelt(i,:), MeltEnergy(i,:), lastpackwater, lastswe,...
                lastsnowdepth] = calcEnergyToMelt(lastswe, lastsnowdepth,...
                packsnowdensity, lastenergy, lastpackwater, LatHeatFreez,...
                WaterDens, SnowMelt(i,:), MeltEnergy(i,:));

            
            % Update water in snowpack
            [Runoff(i,:), lastpackwater] = updatePackWater(a,...
                lastpackwater, lastsnowdepth, S.lw_max, Runoff(i,:),...
                sec_in_ts);
            PackWater(i,:) = lastpackwater;


        %--- Sublimation ---%
        a = lastsnowdepth>0;
        [Sublimation(i,a), Condensation(i,a), lastswe(a),...
            lastsnowdepth(a), lastpackcc(a), packsnowdensity(a),...
            lastpackwater(a)] = calcSublimation(E(a), lastswe(a),...
            lastsnowdepth(a), packsnowdensity(a), lastsnowtemp(a),...
            lastpackcc(a), S.snow_dens_default, Sublimation(i,a),...
            Condensation(i,a), lastpackwater(a), WaterDens);

        
        % update snow:
        b = lastswe > 0;
        lastpackwater(~b) = 0;
        lastalbedo(~b) = S.ground_albedo;
        snowage(~b) = 0;
        lastpacktemp(~b) = 0;
        lastpackcc(~b) = 0;
        lastsnowdepth(~b) = 0;
        lastpacktemp(b) = lastpackcc(b)./ (WaterDens .* Ci .* lastswe(b));
        packsnowdensity(~b) = S.snow_dens_default;


        %--- Apply temperature instability correction ---%
        thres = sec_in_ts/3600 * 0.015; % 15mm for each hour in the time step
        f = find(lastswe < thres & lastswe > 0 & lastpacktemp < (tavg(i,:)));
        lastpacktemp(f) = min(0, tavg(i,f));
        lastpackcc(f) = WaterDens .* Ci .* lastswe(f) .* lastpacktemp(f);

        
        % update outputs:
        PackCC(i,:) = lastpackcc;
        SnowDepth(i,:) = lastsnowdepth;
        SnowWaterEq(i,:) = lastswe;
        b = lastswe>0;
        SnowDensity(i,b) = packsnowdensity(b);

        SW_down(i,:) = Sdn;
        SW_up(i,:) = Sup;
        LW_down(i,:) = Ldn;
        LW_up(i,:) = Lt;
        Q_latent(i,:) = EV;
        Q_sensible(i,:) = H;


    else % If there is not snow:
        lastalbedo = single(ones(1,length(lat))*S.ground_albedo);
        lastswe = single(zeros(1,length(lat)));
        lastsnowdepth = single(zeros(1,length(lat)));
        packsnowdensity = single(ones(1,length(lat))*S.snow_dens_default);
        lastpackwater = single(zeros(1,length(lat)));
        lastpackcc = single(zeros(1,length(lat)));
    end % end if there is snow             

    %disp(cal(i,:))
end % end days


%--- prepare outputs ---%
    SnowWaterEq     = single(SnowWaterEq)*WaterDens;
    SFE             = single(SFE)*WaterDens;
    SnowMelt        = single(SnowMelt)*WaterDens;
    Sublimation     = single(Sublimation)*WaterDens;
    Condensation    = single(Condensation)*WaterDens;
    SnowDepth       = single(SnowDepth)*1000;    
    SnowDensity     = single(SnowDensity);
    Runoff          = single(Runoff)*WaterDens;
    RaininSnow      = single(RaininSnow)*WaterDens;
    RefrozenWater   = single(RefrozenWater)*WaterDens;
    PackWater       = single(PackWater)*WaterDens;
    Albedo          = single(Albedo); %Albedo(SnowYN == 0) = NaN;
    SnowTemp        = single(SnowTemp);
    PackCC          = single(PackCC);
    CCsnowfall      = single(CCsnowfall);
    CCenergy        = single(CCenergy);

    Energy          = single(Energy); 
    MeltEnergy      = single(MeltEnergy); 
    LW_up           = single(LW_up);
    SW_up           = single(SW_up);
    LW_down         = single(LW_down);
    SW_down         = single(SW_down);
    Q_latent        = single(Q_latent);
    Q_sensible      = single(Q_sensible);
    Q_precip        = single(Q_precip);

end        
  