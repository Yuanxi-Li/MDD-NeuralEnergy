function [tauh_CaR_MSN, taum_CaR_MSN] = tau_CaR_MSN(V)
V_CaR_Data = [-30:10:20];  %from data
tauh_CaR_Data = [100 65 35 30 20 20];%from data
if V < V_CaR_Data(1)
    tauh_CaR_MSN = tauh_CaR_Data(1);
elseif V > V_CaR_Data(end)
    tauh_CaR_MSN = tauh_CaR_Data(end);
else
    tauh_CaR_MSN = tauh_CaR_Data(V_CaR_Data==floor(-V/10)*10);
%     tauh_CaR_MSN = interp1(V_CaR_Data, tauh_CaR_Data, V);
end
taum_CaR_MSN = 1.7;
end