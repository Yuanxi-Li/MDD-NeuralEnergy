function [y0] = Initial_PFC_Inter_CB(Vm)
y = zeros(5,1);
y(1) = Vm;

y(4) = 0.002;
[infm_PFC_Inter_CB_Na,...
        tauh_PFC_Inter_CB_Na, infh_PFC_Inter_CB_Na] = inftau_PFC_Inter_CB_mh_Na(y(1));
y(2) = infh_PFC_Inter_CB_Na;

[taun_PFC_Inter_CB_K, infn_PFC_Inter_CB_K] = inftau_PFC_Inter_CB_n_K(y(1));
y(3) = infn_PFC_Inter_CB_K;


[tauh_PFC_Inter_CB_h, infh_PFC_Inter_CB_h] = inftau_PFC_Inter_CB_h_h(y(1));
y(5) = infh_PFC_Inter_CB_h;

y0 = y;

end