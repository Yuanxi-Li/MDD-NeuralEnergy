function [taun_PFC_Inter_CB_K, infn_PFC_Inter_CB_K] = inftau_PFC_Inter_CB_n_K(Vm)
alphan_PFC_Inter_CB_K = -0.01 * (Vm + 34) / (exp(-0.1 * (Vm + 34)) - 1);
betan_PFC_Inter_CB_K = 0.125 * exp(-(Vm + 44) / 80);
[infn_PFC_Inter_CB_K ,taun_PFC_Inter_CB_K] = alpha_beta_solve(alphan_PFC_Inter_CB_K, betan_PFC_Inter_CB_K);
end