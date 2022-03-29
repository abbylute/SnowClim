% Run the SnowClim Model
% ---------------------------

% This script provides an example of how to run the SnowClim Model.
% This script
%  - sets default parameter values
%  - imports example 4-hourly forcing data for water year 2002 for an area 
%    in Central Idaho
%  - runs the SnowClim Model
%  - plots some results
% This script can be modified to work with other parameters, time periods,
% and forcing datasets.


% Start in the SnowClim-Model directory
if ~endsWith(pwd, 'SnowClim-Model')
    error('Please navigate to the SnowClim-Model directory.')
end

addpath('code/')


%--- specify parameters ---%

% to use default values, set parameter equal to []. Parameters hours_in_ts 
% and cal have no default in createParameterFile.m and require the user to 
% specify values.
hours_in_ts = 4;        % number of hours in model run timestep
cal = datevec(datetime(2001,10,1,0,0,0):hours(hours_in_ts):...
    datetime(2002,9,30,23,0,0));
                        % cal specifies the time period to run the model 
                        % for. cal should be 2d array with columns 
                        % corresponding to year, month, day, hour, minute, 
                        % and second and a row for each time step. The 
                        % number of rows in cal should match the number of 
                        % rows in the climate forcing datasets. 
stability = [];         % apply stability correction to turbulent fluxes? 
                        % (0 or 1). Default = 1.
windHt = [];            % height of wind measurements (m). Default = 10.
tempHt = [];            % height of temperature measurements (m). Default =
                        % 2.
snowoff_month = [];     % combined with snowoff_day, specifies the day of 
                        % the year on which all snow should be scraped off 
                        % to avoid making snow towers. (1 to 12). Default =
                        % 9.
snowoff_day = [];       % combined with snowoff_month, specifies the day of
                        % the year on which all snow should be scraped off 
                        % to avoid making snow towers. (1 to 31). Default =
                        % 1.
albedo_option = [];     % albedo algorithm (1=Essery, 2=Tarboton, 3=VIC). 
                        % Default = 2.
max_albedo = [];        % maximum albedo (0 to 1). Default = 0.85.
z_0 = [];               % momentum resistance parameter (m). Default = 
                        % 0.00001.
z_h = [];               % heat and vapor roughness length (m). Default = 
                        % z_0/10.
lw_max = [];            % maximum liquid water content (fraction (0 to 1)  
                        % of snow depth). Default = 0.1.
Tstart = [];            % cold content threshold at which to start energy  
                        % tax (kJ m-2). Default = 0.
Tadd = [];              % cold content range to tax (kJ m-2). Default = 
                        % -10000.
maxtax = [];            % maximum tax to apply to surface energy (0 to 1). 
                        % Default = 0.9.
E0 = [];                % windless heat exchange coefficient (W m-2 K-1). 
                        % Default = 1.
E0_app = [];            % windless heat exchange coefficient flux 
                        % application (1= apply only to sensible heat flux,
                        % 2 = apply to sensible and latent heat fluxes).
                        % Default = 1.
E0_stable = [];         % windless heat exchange coefficient stability 
                        % condition (1 = apply coefficent under any
                        % conditions, 2 = apply only under stable
                        % atmospheric conditions). Default = 2.
Ts_add = [];            % snow surface temperature augmentation (degrees 
                        % C). Default = 2;
smooth_hr = [];         % surface energy flux smoothing window (hours). 
                        % Default = 12.
ground_albedo = [];     % bare ground surface albedo (0 to 1). Default = 
                        % 0.25.
snow_emis = [];         % snow surface emissivity (0 to 1). Default = 0.98.
snow_dens_default = []; % default snowpack density (kg m-3). Default = 250.
G = [];                 % ground heat flux (kJ m-2 s-1). Default = 0.0020.
parameterfilename = 'parameters.mat'; % matfile to save parameters in

createParameterFile(cal, hours_in_ts, stability, windHt,...
    tempHt, snowoff_month, snowoff_day, albedo_option, max_albedo, z_0,...
    z_h, lw_max, Tstart, Tadd, maxtax, E0, E0_app, E0_stable,...
    Ts_add, smooth_hr, ground_albedo, snow_emis, snow_dens_default, G,...
    parameterfilename)


%--- specify forcing data ---%

% lat should be a 1 x n array with n= number of columns in the forcing data
% time varying inputs should be 2d arrays with dimensions time x space
data_dir = 'example_forcings/';

latlonelev = matfile([data_dir,'lat_lon_elev.mat']); 
% latitude in decimal degrees
lat = latlonelev.lat;  
% longitude in decimal degrees
lon = latlonelev.lon;   

% downward longwave radiation at the surface (kJ m-2 hr-1)
lrad = matfile([data_dir,'lrad.mat']);    
lrad = lrad.lrad;

% downward shortwave radiation at the surface (kJ m-2 hr-1)
solar = matfile([data_dir,'solar.mat']);   
solar = solar.solar;   

% average near-surface air temperature (C)
tavg = matfile([data_dir,'tavg.mat']); 
tavg = tavg.tavg; 

% accumulated precipitation (m)
ppt = matfile([data_dir,'ppt.mat']);
ppt = ppt.ppt; 

% wind speed (m s-1)
vs = matfile([data_dir,'vs.mat']);
vs = vs.vs;         

% air pressure (hPa)
psfc = matfile([data_dir,'psfc.mat']);
psfc = psfc.psfc;       

% specific humidity (kg kg-1)
huss = matfile([data_dir,'huss.mat']);
huss = huss.huss;  

% relative humidity (0 to 100)
relhum = matfile([data_dir,'relhum.mat']);
relhum = relhum.relhum;  
% Alternatively, relhum can be calculated as:
% relhum = calcRH(huss, tavg, psfc);
            
% dewpoint temperature (C)
tdmean = matfile([data_dir,'tdmean.mat']);
tdmean = tdmean.tdmean; 
% Alternatively, tdmean can be calculated as:
% tdmean = calcDewpoint(huss,repmat(elev,size(huss,1),1))-273.15; 
% tdmean(tdmean>tavg) = tavg(tdmean>tavg);



%--- run snow model ---%

[SnowMelt, SnowWaterEq, SFE, SnowDepth, SnowDensity,...
    Sublimation, Condensation, SnowTemp, MeltEnergy, Energy, Albedo,...
    SnowYN, RaininSnow, Runoff, RefrozenWater, PackWater, LW_down, LW_up,...
    SW_down, SW_up ,Q_latent, Q_sensible, Q_precip, PackCC, CCenergy,...
    CCsnowfall]...
    = snowclim_model(lat, lrad, tavg, ppt, solar, tdmean,...
    vs, relhum, psfc, huss, parameterfilename);


%--- plot some results ---%

% 1. Plot the snow water equivalent (SWE) niveograph

% calculate spatially average SWE
saswe = mean(SnowWaterEq,2);
% plot
figure(1);clf;
plot(datetime(cal),saswe); hold on;
ylabel('SWE (mm)');
title('Spatially averaged SWE for part of Central Idaho, Water Year 2002');


% 2. Map April 1 Snow Depth

% find april 1:
apr1 = find(cal(:,2)==4 & cal(:,3)==1 & cal(:,4)==0);
% reorganize snow depth
[b,i] = sort(lon);
nlon = length(unique(lon));
depth = SnowDepth(apr1,i);
depth = reshape(depth,[],nlon);
% plot
figure(2);clf;
imagesc(depth);
cb = colorbar();
set(gca,'xtick',(1:10:101));
xt = unique(lon);
set(gca,'xticklabels',xt(1:10:101));
set(gca,'ytick',(1:10:51));
yt = unique(lat);
set(gca,'yticklabels',yt(1:10:51));
title('April 1 Snow Depth, Water Year 2002, Central Idaho');
ylabel(cb,'Snow Depth (mm)');
xlabel('longitude')
ylabel('latitude')



% 3. Plot daily energy fluxes

% convert kJ/m2/timestep to W/m2
% convert to daily averages
% calculate net fluxes
% convert to spatial averages
kj_per_ts_to_w = 1000 /hours_in_ts* 0.00027777777;

% net SW
SW_down(SnowYN==0) = NaN;
SW_downw = SW_down * kj_per_ts_to_w;
SW_downw = to_daily_means(SW_downw, hours_in_ts);

SW_up(SnowYN==0) = NaN;
SW_upw = SW_up * kj_per_ts_to_w;
SW_upw = to_daily_means(SW_upw, hours_in_ts);

SW_net = nanmean(SW_downw,2) - nanmean(SW_upw,2);

% net LW
LW_down(SnowYN==0) = NaN;
LW_downw = LW_down * kj_per_ts_to_w;
LW_downw = to_daily_means(LW_downw, hours_in_ts);

LW_up(SnowYN==0) = NaN;
LW_upw = LW_up * kj_per_ts_to_w;
LW_upw = to_daily_means(LW_upw, hours_in_ts);

LW_net = nanmean(LW_downw,2) - nanmean(LW_upw,2);

% sensible 
Q_sensible(SnowYN==0) = NaN;
Q_sensiblew = Q_sensible * kj_per_ts_to_w;
Q_sensiblew = to_daily_means(Q_sensiblew, hours_in_ts);
Q_sensiblew = nanmean(Q_sensiblew,2);

% latent
Q_latent(SnowYN==0) = NaN;
Q_latentw = Q_latent * kj_per_ts_to_w;
Q_latentw = to_daily_means(Q_latentw, hours_in_ts);
Q_latentw = nanmean(Q_latentw,2);

% precip
Q_precip = Q_precip + CCsnowfall; % combine liquid and solid heat fluxes
Q_precip(SnowYN==0) = NaN;
Q_precipw = Q_precip * kj_per_ts_to_w;
Q_precipw = to_daily_means(Q_precipw, hours_in_ts);
Q_precipw = nanmean(Q_precipw,2);

% ground
m = matfile(parameterfilename);
m = m.S;
G = m.G; % ground heat flux parameter in kJ m-2 s-1
G = ones(size(Energy)) * G; % create a timeseries
G(SnowYN==0) = NaN;
G = G * 60 * 60 * hours_in_ts; % kJ m-2 timestep-1
Q_groundw = G * kj_per_ts_to_w;
Q_groundw = to_daily_means(Q_groundw, hours_in_ts);
Q_groundw = nanmean(Q_groundw,2);

% net energy
Energy(SnowYN==0) = NaN;
Q_energyw = Energy * kj_per_ts_to_w;
Q_energyw = to_daily_means(Q_energyw, hours_in_ts);
Q_energyw = nanmean(Q_energyw,2);

% create daily calendar to plot with
dcal = unique(datetime(cal(:,1:3)));

% plot the daily energy fluxes, spatially averaged over locations with snow
% cover
figure(3);clf;
plot(dcal, SW_net); hold on;
plot(dcal, LW_net);
plot(dcal, Q_sensiblew);
plot(dcal, Q_latentw);
plot(dcal, Q_precipw);
plot(dcal, Q_groundw);
plot(dcal, Q_energyw,'k');
set(gca,'XLim',[dcal(1),dcal(290)]);
ylabel('Energy Flux (W/m^{2})');
legend('Net SW','Net LW','Sensible','Latent','Precip','Ground','Net Energy', 'Location','NorthWest')
title('Water Year 2002 Daily Energy Fluxes')

