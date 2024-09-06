function [tauq_PFC_Pyra_Proximal_KS, infq_PFC_Pyra_Proximal_KS,...
    taur_PFC_Pyra_Proximal_KS, infr_PFC_Pyra_Proximal_KS] = inftau_PFC_Pyra_Proximal_qr_KS(Vm)
infq_PFC_Pyra_Proximal_KS = 1/(1+exp(-(Vm+34)/6.5));
tauq_PFC_Pyra_Proximal_KS = 8/(exp(-((Vm+55)/30))+exp((Vm+55)/30));
infr_PFC_Pyra_Proximal_KS = 1/(1+exp((Vm+65)/6.6));
taur_PFC_Pyra_Proximal_KS = 100/(1+exp(-((Vm+65)/6.8)))+100;

end