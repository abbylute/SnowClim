function [incangle] = calcInclinationAngle(lat, month, day)
% calculate the solar inclination angle as a function of latitude and 
% time of year

% if day is not provided, use midday of month
if nargin < 3
    day = 15;
end

doy = datenum(2013,month,day) - datenum(2013,1,0);
decl_angle = 23.45 * sin(2*pi*(284+doy)/365);
incangle = 90 - abs(lat - decl_angle);