function [tauh_NaF_MSN, taum_NaF_MSN] = tau_NaF_MSN(V);
V_NaF_Data = [-100:10:50];  %from data
tauh_NaF_Data = [1.3 1.3 1.3 1.3 1.3 1.3 1.3 1.3 0.85 0.5 0.45 0.32 0.30 0.28 0.28 0.28];%from data
taum_NaF_Data = [0.06 0.06 0.07 0.09 0.11 0.13 0.20 0.32 0.16 0.15 0.12 0.08 0.06 0.06 0.06 0.06]; %from data
if V < V_NaF_Data(1)
    tauh_NaF_MSN = tauh_NaF_Data(1);
    taum_NaF_MSN = taum_NaF_Data(1);
elseif V > V_NaF_Data(end)
    tauh_NaF_MSN = tauh_NaF_Data(end);
    taum_NaF_MSN = taum_NaF_Data(end);
else
    tauh_NaF_MSN = tauh_NaF_Data(V_NaF_Data==floor(V/10)*10);
    taum_NaF_MSN = taum_NaF_Data(V_NaF_Data==floor(V/10)*10);
%     tauh_NaF_MSN = interp1(V_NaF_Data, tauh_NaF_Data, V);
%     taum_NaF_MSN = interp1(V_NaF_Data, taum_NaF_Data, V);
end
end