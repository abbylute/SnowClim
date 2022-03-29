function [H, E, EV] = calcTurbulentFluxes(stability, ...
    windHt, z_0, tempHt, z_h, k, vs, lastsnowtemp, tavg, ...
    Ca, psfc, huss, E0_value, E0_app, E0_stable, sec_in_ts)

      
 % Essery et al., (2013), option 1: Richardson number parameterization
        
        % Set parameters:
        R = 287;           % gas constant for dry air (J K-1 kg-1)
        g = 9.80616;       % acceleration due to gravity (m s-2)
        c = 5;             % from Essery et al., (2013) table 6, taken from Louis (1979)
        E0 = E0_value/1000; % windless exchange coefficient (approx 2 W/m2/K (converted to kJ/m2/K/s), based on Brown et al., 2006, drawing on the work of Jordan, 1991 and 1999)
        
        % convert air pressure from hPa to Pa:
        psfcpa = psfc*100;
        
        % calculate air density (kg/m3):
        pa = psfcpa ./ (R .* (tavg + 273.15));
        
        % calculate vapor densities (kg/m3):
        rhoa = huss;
        rhos = calcSpecificHumidity(lastsnowtemp,psfc);
        
        % calculate exchange coefficient for neutral conditions CHN (i.e. for Rib = 0)
        CHN = k^2 .* (log(windHt/z_0))^(-1) .* (log(tempHt/z_h))^(-1);

        if stability
            % calculate the bulk Richardson number:
            Rib = (g * windHt * (tavg - lastsnowtemp))./ ((tavg + 273.15) .* vs.^2);

            % calculate FH as a function of Rib:
            FH = ones(size(tavg))*NaN;
            % for the unstable case:
            b = Rib < 0; 
            FH(b) = 1 - ((3*c*Rib(b)) ./ (1 + 3*c^2*CHN* (-Rib(b)*windHt/z_0).^.5));
            % for neutral case:
            b = Rib == 0;
            FH(b) = 1;
            % for the stable case:
            b = Rib > 0;
            FH(b) = (1 + ((2*c*Rib(b)) ./ (1 + Rib(b)).^.5)).^(-1);

            % calculate the exchange coefficient CH:
            CH = FH .* CHN;
        else
            CH = CHN;
        end

        LatHeatVap = calcLatHeatVap(lastsnowtemp); % (kJ/kg)
        LatHeatSub = calcLatHeatSub(lastsnowtemp); % (kJ/kg)

        % calculate windless exchange coefficient
        Ex = zeros(size(vs));
        if E0_stable==1
            Ex = E0;
        elseif E0_stable==2
            Ex(Rib>0) = E0;
        end
        H = -(pa .* Ca .* CH .* vs + Ex) .* (lastsnowtemp - tavg);
        
        if E0_app==1
            E = -(pa .* CH .* vs) .* (rhos - rhoa); % mass flux kg/m2/s
        elseif E0_app==2
            E = -(pa .* CH .* vs + Ex) .* (rhos - rhoa); % mass flux kg/m2/s
        end
        Evap = E .* LatHeatVap;
        Esub = E .* LatHeatSub;
        
        EV = zeros(size(E));
        EV(lastsnowtemp >= 0) = Evap(lastsnowtemp >=0);
        EV(lastsnowtemp < 0) = Esub(lastsnowtemp < 0);
        
        % convert from /second to /timestep
        H = H * sec_in_ts; % kJ/m2/s
        E = E * sec_in_ts; % kg/m2/s
        EV = EV * sec_in_ts; % kJ/m2/s

        % avoid NaNs
        H(vs==0)=0;
        E(vs==0)=0;
        EV(vs==0)=0;
end

