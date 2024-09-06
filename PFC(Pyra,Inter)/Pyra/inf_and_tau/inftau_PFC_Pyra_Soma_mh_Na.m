function [infm_PFC_Pyra_Soma_Na,...
    tauh_PFC_Pyra_Soma_Na, infh_PFC_Pyra_Soma_Na] = inftau_PFC_Pyra_Soma_mh_Na(Vm)
alpham_PFC_Pyra_Soma_Na = -0.1*(Vm+31)/(exp(-0.1*(Vm+31))-1);
betam_PFC_Pyra_Soma_Na = 4*exp(-(Vm+56)/18);
alphah_PFC_Pyra_Soma_Na = 0.07*exp(-(Vm+47)/20);
betah_PFC_Pyra_Soma_Na = 1/(exp(-0.1*(Vm+17))+1);
infm_PFC_Pyra_Soma_Na = alpham_PFC_Pyra_Soma_Na/(alpham_PFC_Pyra_Soma_Na+betam_PFC_Pyra_Soma_Na);
[infh_PFC_Pyra_Soma_Na,tauh_PFC_Pyra_Soma_Na] = alpha_beta_solve(alphah_PFC_Pyra_Soma_Na, betah_PFC_Pyra_Soma_Na);
end