function [tauh_CaN_MSN, taum_CaN_MSN] = tau_CaN_MSN(V);
malpha_CaN_MSN = 0.1157*(V+17.19)/(exp((V+17.19)/15.22)-1);
mbeta_CaN_MSN = 1.15*exp(V/23.82);
taum_CaN_MSN = 1/(malpha_CaN_MSN+mbeta_CaN_MSN);
tauh_CaN_MSN = 23.33;
end