function [tauh_KAf_MSN, taum_KAf_MSN] = tau_KAf_MSN(V);
V_KAf_Data = [-40:10:60];  %from data
taum_KAf_Data = [1.8 1.1 1.0 1.0 0.9 0.8 0.9 0.9 0.9 0.8 0.8]; %from data
if V < V_KAf_Data(1)
    taum_KAf_MSN = taum_KAf_Data(1);
elseif V > V_KAf_Data(end)
    taum_KAf_MSN = taum_KAf_Data(end);
else
    taum_KAf_MSN = taum_KAf_Data(V_KAf_Data==floor(V/10)*10);
%     taum_KAf_MSN = interp1(V_KAf_Data, taum_KAf_Data, V);
end
tauh_KAf_MSN = 4.67;
end