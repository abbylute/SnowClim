function [lw] = calcLongwave(emissivity, temp, lwdown, sec_in_ts)
% Calculate terrestrial upward longwave radiation

    SBconstant = 5.67E-11;%         #  [kJ m-2 K-4 s-1]
    tempK = temp + 273.15;%         #  [degrees K]
    lw = (emissivity.*SBconstant.*tempK.^4) + ((1-emissivity)*(lwdown/sec_in_ts));
    lw = lw .* sec_in_ts;