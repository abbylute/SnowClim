
function relhum=calcRH(huss, temp, pres)
% based on:
% https://earthscience.stackexchange.com/questions/2360/how-do-i-convert-specific-humidity-to-relative-humidity

% might consider using Dai, 2006 instead

% huss: specific humidity (kg/kg)
% temp: air temperature (C)
% pres: surface pressure (mb)


es = 6.112*exp((17.67.*temp)./(temp + 243.5));
e = huss.*pres./(0.378.*huss+0.622);
relhum = 100*e./es;
relhum(relhum>100) = 100;
relhum(relhum<0) = 0;