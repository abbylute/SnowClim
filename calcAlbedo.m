function [A, snowage] = calcAlbedo(albedo_option, groundAlbedo, maxAlbedo,... 
    lastalbedo,newsnowdepth, lastsnowdepth, newswe, lastswe, lastsnowtemp,...
    WaterDens, lat, month, day, snowage, lastpackcc, sec_in_ts)
  
switch albedo_option
    case 1 % Essery et al., (2013) option 1
        
        A = calc_albedo_Essery_opt1(lastalbedo, lastsnowtemp, newswe,...
            lastswe, maxAlbedo, groundAlbedo, WaterDens, sec_in_ts);
     
    case 2 % Utah Snow Model (Tarboton)
        
        [A, snowage]=calc_albedo_tarboton(lastsnowtemp+273.15, snowage,...
            newsnowdepth, lat, month, day, lastswe, groundAlbedo,...
            maxAlbedo,sec_in_ts);
        
    case 3 %VIC

        [A, snowage]=calc_albedo_vic(newsnowdepth, snowage, groundAlbedo,...
            lastsnowdepth, lastalbedo, maxAlbedo, sec_in_ts, lastpackcc);

% if total snow depth is low <0.1 m, then adjust albedo
    z = lastsnowdepth;
    f = z < 0.1;
    r = (1 - z(f)/.1).*exp(-z(f)./0.2);
    r = min(r, 1);
    A(f) = r.*groundAlbedo + (1 - r).*A(f);


A = real(A);
 
end



function A = calc_albedo_Essery_opt1(lastalbedo, lastsnowtemp, newswe,...
        lastswe, maxAlbedo, groundAlbedo, WaterDens, sec_in_ts)

% calculate snow albedo
% based on Essery et al., (2013) option 1 (eqns 26-28).
% Parameter values are from Table 5 in Essery et al., (2013), which are
% based on Douville et al., (1995).

% parameters:
minAlbedo = 0.5;    % minimum snow albedo
So = 10;            % kg/m2 ('critical SWE')
Ta = 10^7;          % s
Tm = 3.6E5;         % s

dt = sec_in_ts;      % time step in seconds
Sf = newswe.*WaterDens/dt;      % snow fall rate (kg m-2 s-1)


% if there's no snow on the ground:
    b = (lastswe == 0);
    A(b) = groundAlbedo;

% if there's fresh snow:
    b = newswe > 0;
    A(b) = lastalbedo(b) + (maxAlbedo - lastalbedo(b)) .* ((Sf(b).*dt)/So);

% for cold snow:
    b = lastswe > 0 & newswe == 0 & lastsnowtemp < -.5;
    A(b) = lastalbedo(b) - dt/Ta; % less change for colder snow

% for melting snow:
    b = lastswe > 0 & newswe == 0 & lastsnowtemp >= -.5;
    A(b) = (lastalbedo(b) - minAlbedo) .* exp(-dt/Tm) + minAlbedo;
    
    
A(A < minAlbedo) = minAlbedo;
A(A > maxAlbedo) = maxAlbedo;






 function [A, snowage]=calc_albedo_vic(newsnowdepth, snowage, groundAlbedo, ...
         lastsnowdepth, lastalbedo, maxAlbedo, sec_in_ts, lastpackcc)
     
 % adapted from the VIC code snow_utility.c

 SNOW_NEW_SNOW_ALB=maxAlbedo; 
 SNOW_ALB_ACCUM_A=0.94;
 SNOW_ALB_ACCUM_B=0.58;
 SNOW_ALB_THAW_A=0.82;
 SNOW_ALB_THAW_B=0.46;
 sec_per_day=24*60*60;
 
 % default the new snow albedo to lastalbedo
 A = lastalbedo;
 
 % new snow case:
 b = newsnowdepth > 0.01 & lastpackcc < 0;
 snowage(b) = 0;
 A(b) = SNOW_NEW_SNOW_ALB;
 
 % aged snow case:
 b = ~b & (lastsnowdepth > 0);
 snowage(b) = snowage(b) + sec_in_ts;
    % if cold:
    c = b & lastpackcc < 0;
    A(c) = SNOW_NEW_SNOW_ALB .* SNOW_ALB_ACCUM_A.^((snowage(c)./sec_per_day).^SNOW_ALB_ACCUM_B);
 
    % if melting:
    c = b & lastpackcc == 0;
    A(c) = SNOW_NEW_SNOW_ALB .* SNOW_ALB_THAW_A.^((snowage(c)./sec_per_day).^SNOW_ALB_THAW_B);
 
 b = lastsnowdepth == 0;
 snowage(b) = 0;
 A(b) = groundAlbedo;


 
 
function [albedo, snowage]=calc_albedo_tarboton(lastsnowtemp, snowage,...
        newsnowdepth, lat, month, day, lastswe, groundAlbedo, maxAlbedo,...
        sec_in_ts)
% adapted from UEB snow model:
% fs.fed.us/rm/boise/publications/watershed/rmrs_1996_tarbotond001.pdf
    
    Cv  = 0.2;
    Cir = 0.5;

    albedo_iro = 0.65;
    albedo_vo  = 2*maxAlbedo-albedo_iro;
    if albedo_vo>1
        d = albedo_vo-1;
        albedo_iro=albedo_iro+d;
        albedo_vo=1;
    end

    % for new snow depths > 0.01 m
    b = newsnowdepth > 0.01;
        albedo(b) = albedo_vo/2 + albedo_iro/2;
        snowage(b) = 0;
        inc = calcInclinationAngle(lat(b), month, day) * pi/180;
        c = cos(inc) < .5;
            fpsi = 0.5 * (3./(1+4*cos(inc(c))) - 1);
            extra = zeros(size(inc));
            extra(c) = 0.2 * fpsi * (1-albedo_vo) + 0.2*fpsi*(1-albedo_iro);
            albedo(b) = albedo(b) + extra;

    % for new snow depths <= 0.01 m albedo is a function of snow age
    b = ~b;
        r1 = exp(5000 * (1/273.16 - 1./lastsnowtemp(b)));
        r2 = min(r1.^10, 1);
        r3 = 0.03;

        inc_age = (r1+r2+r3)/1e6 * sec_in_ts;

        snowage(b) = snowage(b) + inc_age;

        Fage = snowage(b) ./ (1+snowage(b));

        albedo_1 = (1 - Cv.*Fage) .* albedo_vo;
        albedo_2 = (1 - Cir.*Fage) .* albedo_iro;

        albedo(b) = albedo_1/2 + albedo_2/2;
        inc = calcInclinationAngle(lat(b), month, day) * pi/180;
        c = cos(inc) < .5;
            fpsi = 0.5 * (3./(1+4*cos(inc(c)))-1);
            extra = zeros(size(inc));
            extra(c) = 0.2 .* fpsi .* (1-albedo_1(c)) + 0.2 .* fpsi .* (1-albedo_2(c));
            albedo(b) = albedo(b) + extra;

    % for no snow case
    b = lastswe == 0;
    albedo(b) = groundAlbedo;
    snowage(b) = 0;

