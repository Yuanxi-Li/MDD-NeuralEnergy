function [tauh_PFC_Inter_CB_h, infh_PFC_Inter_CB_h] = inftau_PFC_Inter_CB_h_h(Vm)
infh_PFC_Inter_CB_h = 1 / (1 + exp((Vm + 80)/10));
tauh_PFC_Inter_CB_h = 200 / (exp((Vm + 70)/20) + exp(-(Vm+70)/20)) + 5;
end