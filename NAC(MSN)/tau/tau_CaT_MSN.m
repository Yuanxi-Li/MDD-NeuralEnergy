function [tauh_CaT_MSN, taum_CaT_MSN] = tau_CaT_MSN(V);
V_CaT_Data = [-65:5:10];  %from data
tauh_CaT_Data = [382 208 162 129 119 107 107 107 108 109 109 110 110 110 110 110];%from data
taum_CaT_Data = [20.2 20.2 13.1 8.7 6.8 5.6 4.4 3.8 3.6 3.3 3.6 3.6 3.3 3.3 3.3 3.3]; %from data
if V < V_CaT_Data(1)
    tauh_CaT_MSN = tauh_CaT_Data(1);
    taum_CaT_MSN = taum_CaT_Data(1);
elseif V > V_CaT_Data(end)
    tauh_CaT_MSN = tauh_CaT_Data(end);
    taum_CaT_MSN = taum_CaT_Data(end);
else
    tauh_CaT_MSN = tauh_CaT_Data(V_CaT_Data==floor(V/5)*5);
    taum_CaT_MSN = taum_CaT_Data(V_CaT_Data==floor(V/5)*5);
%     tauh_CaT_MSN = interp1(V_CaT_Data, tauh_CaT_Data, V);
%     taum_CaT_MSN = interp1(V_CaT_Data, taum_CaT_Data, V);
end
end