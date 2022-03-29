function [snow]=calcPhase(temp, precip, rh, binary)

% this function is based from the bivariate logistic regression model 
% presented in Jennings et al., 2018.  It uses air temperature and relative
% humidity to estimate precipitation phase. When binary=True, precipitation
% is classified as snow if the probability is >=50%, and rain otherwise. 
% Coefficient values (a, b, g) are drawn from the supplementary material in 
% Jennings et al., 2018, table 2.

% temp should be in K

a = -10.04;
b = 1.41;
g = 0.09;

temp=temp-273.15;

psnow = 1./(1 + exp(a + b*temp + g*rh));

if binary == true
    psnow(psnow<.5) = 0;
    psnow(psnow>=.5) = 1;
end
snow=single(psnow.*precip);
