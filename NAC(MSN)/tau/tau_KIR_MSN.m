function taum_KIR_MSN = tau_KIR_MSN(V);
V_KIR_Data = [-100:10:50];  %from data
taum_KIR_Data = [3.7313 4.0000 4.7170 5.3763 6.0606 6.8966 7.6923 7.1429 5.8824 4.4444 4.0000 4.0000 4.0000 4.0000 4.0000 4.0000]; %from data
if V < V_KIR_Data(1)
    taum_KIR_MSN = taum_KIR_Data(1);
elseif V > V_KIR_Data(end)
    taum_KIR_MSN = taum_KIR_Data(end);
else
    taum_KIR_MSN = taum_KIR_Data(V_KIR_Data==floor(V/10)*10);
%     taum_KIR_MSN = interp1(V_KIR_Data, taum_KIR_Data, V);
end
end