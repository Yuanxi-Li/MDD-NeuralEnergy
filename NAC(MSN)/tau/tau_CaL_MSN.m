function [tauh_CaL_MSN, taum_CaL_MSN] = tau_CaL_MSN(V);
malpha_CaL_MSN = 0.1194*(V+8.124)/(exp((V+8.124)/9.005)-1);
mbeta_CaL_MSN = 2.97*exp(V/31.4);
taum_CaL_MSN = 1/(malpha_CaL_MSN+mbeta_CaL_MSN);
tauh_CaL_MSN = 14.77;
end