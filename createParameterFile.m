function[] = createParameterFile(cal, hours_in_ts, stability, windHt,...
    tempHt, snowoff_month, snowoff_day, albedo_option, max_albedo, z_0,...
    z_h, lw_max, Tstart, Tadd, maxtax, E0, E0_app, E0_stable, Ts_add,...
    smooth_hr, ground_albedo, snow_emis, snow_dens_default, G,...
    parameterfilename)

% write parameters to structure for snow model to read

% default values are those used in application of the snow model to the
% western United States, some of which were determined through calibration 
% at SNOTEL stations (Lute et al., in prep).

S = struct();
if isempty(cal)
    error('Argument cal is required.')
end
S.cal = cal; %datevec(datetime(2000,10,1,0,0,0):hours(S.hours_in_ts):datetime(2013,9,30,23,0,0));

S.hours_in_ts = hours_in_ts;

if isempty(stability)
    stability = 1;
end
S.stability = stability;

if isempty(windHt)
    windHt = 10;
end
S.windHt = windHt;

if isempty(tempHt)
    tempHt = 2;
end
S.tempHt = tempHt;

if isempty(snowoff_month)
    snowoff_month = 9;
end
S.snowoff_month = snowoff_month;

if isempty(snowoff_day)
    snowoff_day = 1;
end
S.snowoff_day = snowoff_day;

if isempty(albedo_option)
    albedo_option = 2;
end
S.albedo_option = albedo_option;

if isempty(max_albedo)
    max_albedo = 0.85;
end
S.max_albedo = max_albedo;

if isempty(z_0)
    z_0 = 0.00001;
end
S.z_0 = z_0;

if isempty(z_h)
    z_h = z_0/10;
end
S.z_h = z_h;

if isempty(lw_max)
    lw_max = 0.1;
end
S.lw_max = lw_max;

if isempty(Tstart)
    Tstart = 0;
end
S.Tstart = Tstart;

if isempty(Tadd)
    Tadd = -10000;
end
S.Tadd = Tadd;

if isempty(maxtax)
    maxtax = 0.9;
end
S.maxtax = maxtax;

if isempty(E0)
    E0 = 1;
end
S.E0_value = E0;

if isempty(E0_app)
    E0_app = 1;
end
S.E0_app = E0_app;

if isempty(E0_stable)
    E0_stable = 2;
end
S.E0_stable = E0_stable;

if isempty(Ts_add)
    Ts_add = 2;
end
S.Ts_add = Ts_add;

if isempty(smooth_hr)
    smooth_hr = 12;
end
S.smooth_hr = smooth_hr;

if isempty(ground_albedo)
    ground_albedo = 0.25;
end
S.ground_albedo = ground_albedo;

if isempty(snow_emis)
    snow_emis = 0.98;   % snow emissivity (from Snow and Climate, Armstrong + Brun eds, pg 58)
end
S.snow_emis = snow_emis;

if isempty(snow_dens_default)
    snow_dens_default = 250;  % default snow density (kg/m3) (from Essery et al., (2013), snow compaction option 2, based on Cox et al., (1999))
end
S.snow_dens_default = snow_dens_default;   

if isempty(G)
    G = 173/86400;  % ground conduction (kJ/m2/s), from Walter et al., (2005)
end
S.G = G;

save(parameterfilename,'S');

