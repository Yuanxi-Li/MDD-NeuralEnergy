function [tauh_NaP_MSN, taum_NaP_MSN] = tau_NaP_MSN(V);
V_NaP_Data = -100:10:40;  %from data
tauh_NaP_Data = [4500 4750 5200 6100 6300 5000 4250 3500 3000 2700 2500 2100 2100 2100 2100]; %from data
if V < V_NaP_Data(1)
    tauh_NaP_MSN = tauh_NaP_Data(1);
elseif V > V_NaP_Data(end)
    tauh_NaP_MSN = tauh_NaP_Data(end);
else
    tauh_NaP_MSN = tauh_NaP_Data(V_NaP_Data==floor(V/10)*10);
%     tauh_NaP_MSN = interp1(V_NaP_Data, tauh_NaP_Data, V);
end
% taum
if V < -40
    taum_NaP_MSN = 0.025+0.14*exp((V+40)/10);
elseif V >= -40
    taum_NaP_MSN = 0.02+0.145*exp(-(V+40)/10);
end

