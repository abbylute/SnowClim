function [sh] = calcSpecificHumidity(Td, P)
% Td = dewpoint temperature (C)
% P = surface pressure (mb)

e = 6.112.*exp((17.67.*Td)./(Td+243.5));
sh = (0.622.*e)./(P - (0.378.*e));

% sources:
% https://archive.eol.ucar.edu/projects/ceop/dm/documents/refdata_report/eqns.html
% the above cites Bolton (1980):
% https://journals.ametsoc.org/mwr/article/108/7/1046/62205/The-Computation-of-Equivalent-Potential
end