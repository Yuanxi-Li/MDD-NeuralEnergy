function [tauh_KRP_MSN, taum_KRP_MSN] = tau_KRP_MSN(V);
V_KRP_Data = [-100:5:50];  %from data
tauh_KRP_Data = [7000.0 7000.0 7000.0 7000.0 7000.0 7000.0 7000.0 7000.0 7000.0 7000.0 7000.0 7000.0 ...
     7000.0 7000.0 7000.0 7000.0 7000.0 7000.0 7000.0 7000.0 6742.5 6000.0 4740.2 3500.0 2783.3 2500.0 ...
     2336.3 2200.0 2083.5 2000.0 2000.0];%from data
taum_KRP_Data = [40 45 48.8 55 64.4 75 83.9 90 93.5 95 95.4 97 99.2 95 79.7 60 44.5 ...
    35 29.3 25 20 15 11.6 10 9.6 10 10.5 10 8 5 5]; %from data
if V < V_KRP_Data(1)
    tauh_KRP_MSN = tauh_KRP_Data(1);
    taum_KRP_MSN = taum_KRP_Data(1);
elseif V > V_KRP_Data(end)
    tauh_KRP_MSN = tauh_KRP_Data(end);
    taum_KRP_MSN = taum_KRP_Data(end);
else
    tauh_KRP_MSN = tauh_KRP_Data(V_KRP_Data==floor(V/5)*5);
    taum_KRP_MSN = taum_KRP_Data(V_KRP_Data==floor(V/5)*5);
    %
%     tauh_KRP_MSN = interp1(V_KRP_Data, tauh_KRP_Data, V);
%     taum_KRP_MSN = interp1(V_KRP_Data, taum_KRP_Data, V);
end
end