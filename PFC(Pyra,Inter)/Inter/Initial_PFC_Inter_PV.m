function [y0] = Initial_PFC_Inter_PV(V)
y = zeros(3,1);
y(1) = V;
[infm_PFC_Inter_PV_Na,...
        tauh_PFC_Inter_PV_Na, infh_PFC_Inter_PV_Na] = inftau_PFC_Inter_PV_mh_Na(y(1));
y(2) = infh_PFC_Inter_PV_Na;
[taun_PFC_Inter_PV_K, infn_PFC_Inter_PV_K] = inftau_PFC_Inter_PV_n_K(y(1));
y(3) = infn_PFC_Inter_PV_K;
y0 = y;
end