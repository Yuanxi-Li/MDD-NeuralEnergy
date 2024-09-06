function [infm_PFC_Inter_PV_Na,...
    tauh_PFC_Inter_PV_Na, infh_PFC_Inter_PV_Na] = inftau_PFC_Inter_PV_mh_Na(Vm)
alpham_PFC_Inter_PV_Na = -0.1 * (Vm + 35) / (exp(-0.1 * (Vm + 35)) - 1);
betam_PFC_Inter_PV_Na = 4 * exp(-(Vm + 60) / 18);
alphah_PFC_Inter_PV_Na = 0.07 * exp(-(Vm + 58) / 20);
betah_PFC_Inter_PV_Na = 1 / (exp(-0.1 * (Vm + 28)) + 1);
infm_PFC_Inter_PV_Na = alpham_PFC_Inter_PV_Na / (alpham_PFC_Inter_PV_Na + betam_PFC_Inter_PV_Na);
[infh_PFC_Inter_PV_Na ,tauh_PFC_Inter_PV_Na] = alpha_beta_solve(alphah_PFC_Inter_PV_Na, betah_PFC_Inter_PV_Na);

end