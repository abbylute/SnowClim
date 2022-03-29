function [Td]=calcDewpoint(shv,Z);

pres = 1013.25*((293-0.0065*Z)/293).^5.26; % CIMIS
shv=shv.*pres/.622;
e1=log(shv/6.112);
Td=243.5*e1./(17.67-e1);
Td=Td+273.15; %if in Kelvin
