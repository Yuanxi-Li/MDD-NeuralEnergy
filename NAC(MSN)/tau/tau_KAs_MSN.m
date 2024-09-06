function [tauh_KAs_MSN, taum_KAs_MSN] = tau_KAs_MSN(V);
taum_KAs_MSN = 0.378+9.91*exp(-((V+34.3)/30.1)^2);
halpha_KAs_Soma = exp(-(V+90.96)/29.01);
hbeta_KAs_Soma = exp((V+90.96)/100);
tauh_KAs_MSN = 1097.4/(halpha_KAs_Soma+hbeta_KAs_Soma);
end