function[runoff, lastpackwater] = updatePackWater(a, lastpackwater,...
    lastsnowdepth, lw_max, runoff, sec_in_ts)

    % Update snowpack liquid water content 

    % irreducible water (1% from Marsh, 1991, referenced in Kinar and Pomeroy 2015)
    min_water = .01 .*lastsnowdepth; 
    % maximum water content (10% as referenced in Kinar and Pomeroy 2015)
    % this could be higher, maybe .15 or even .2 (from snow principles pg69)
    max_water = lw_max .*lastsnowdepth; 
    
    % if packwater <= min_water, do nothing

    % if packwater > min_water & < max_water, allow gravity drainage at
    % a rate of 10 cm/hr
    b = a & lastpackwater > min_water & lastpackwater <= max_water;
        % drainage rate of 10 cm/hr, converted to m/s
        waterrate = ones(size(lastpackwater(b))).*2.7778e-05; 
        % translate to volume (m):
        % divide sec_in_ts by 2 since we apply updatePackWater twice each time step
        waterdrainage = waterrate.*(sec_in_ts/2);
        % don't let drainage take packwater below minimum
        waterdrainage = min(waterdrainage, lastpackwater(b)-min_water(b));
        runoff(:,b) = runoff(:,b) + waterdrainage;
        lastpackwater(b) = lastpackwater(b) - waterdrainage;
        
    % if packwater > max_water, drain all excess
    b = a & lastpackwater > max_water;
    runoff(:,b) = runoff(:,b) + (lastpackwater(b) - max_water(b));
    lastpackwater(b) = max_water(b);
    
end
