function [taun_PFC_Pyra_Soma_K, infn_PFC_Pyra_Soma_K] = inftau_PFC_Pyra_Soma_n_K(Vm)
alphan_PFC_Pyra_Soma_K = -0.01*(Vm+34)/(exp(-0.1*(Vm+34))-1);
betan_PFC_Pyra_Soma_K = 0.125*exp(-(Vm+44)/80);
[infn_PFC_Pyra_Soma_K,taun_PFC_Pyra_Soma_K] = alpha_beta_solve(alphan_PFC_Pyra_Soma_K, betan_PFC_Pyra_Soma_K);
end