function[lastpackcc, lastenergy, CCenergy] = calcEnergyToCC(lastpackcc,...
    lastenergy, CCenergy)
    % distribute excess energy to snowpack cold content

    % 1. if cold content = 0 and lastenergy > 0, then nothing changes

    % 2. if -lastpackcc > 0 but less than or equal to lastenergy, then 
    b = -lastpackcc > 0 & -lastpackcc <= lastenergy;
    lastenergy(b) = lastenergy(b) + lastpackcc(b);
    lastpackcc(b) = 0;
    CCenergy(:,b) = min(0,lastenergy(b));
    % and the rest goes to melt

    % 3. if -lastpackcc > lastenergy, then energy warms the snow, and
    % some cc can go to refreezing.
    % includes cases with cc=0 and lastenergy<0, also any case with 
    % lastenergy<0
    b = -lastpackcc > lastenergy; 
    lastpackcc(b) = lastpackcc(b) + lastenergy(b);
    CCenergy(:,b) = min(0,lastenergy(b));
    lastenergy(b) = 0;
    
end