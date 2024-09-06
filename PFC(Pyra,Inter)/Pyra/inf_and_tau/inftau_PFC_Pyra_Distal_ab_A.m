function [taua_PFC_Pyra_Distal_A, infa_PFC_Pyra_Distal_A,...
    taub_PFC_Pyra_Distal_A, infb_PFC_Pyra_Distal_A] = inftau_PFC_Pyra_Distal_ab_A(Vm)


infa_PFC_Pyra_Distal_A = 1/(1+exp(-(Vm+60)/8.5));
taua_PFC_Pyra_Distal_A = 0.37+1/(exp((Vm+35.8)/19.7)+exp(-(Vm+79.7)/12.7));
infb_PFC_Pyra_Distal_A = 1/(1+exp((Vm+78)/6));
taub_PFC_Pyra_Distal_A = 19+1/(exp((Vm+46)/5)+exp((Vm+238)/(-37.5)));

end