function [y0] = Initial_PFC_Pyra(V_Soma,V_Proximal,V_Distal, Dopamine_Ratio);
y = zeros(12,1);
y(1) = V_Soma;
y(2) = V_Proximal;
y(3) = V_Distal;
y(7) = 0.000667;
[infm_PFC_Pyra_Soma_Na,...
    tauh_PFC_Pyra_Soma_Na, infh_PFC_Pyra_Soma_Na] = inftau_PFC_Pyra_Soma_mh_Na(y(1));

y(4) = infh_PFC_Pyra_Soma_Na;
[infn_PFC_Pyra_Soma_K,taun_PFC_Pyra_Soma_K] = inftau_PFC_Pyra_Soma_n_K(y(1));
y(5) = infn_PFC_Pyra_Soma_K;
alpham_PFC_Pyra_Soma_CaN = 0.0056;
betam_PFC_Pyra_Soma_CaN = 0.002;
infm_PFC_Pyra_Soma_CaN = alpham_PFC_Pyra_Soma_CaN*(y(7)^2)/(alpham_PFC_Pyra_Soma_CaN*(y(7)^2)+betam_PFC_Pyra_Soma_CaN);

y(6) = infm_PFC_Pyra_Soma_CaN;

[infm_PFC_Pyra_Proximal_Na,...
    tauh_PFC_Pyra_Proximal_Na, infh_PFC_Pyra_Proximal_Na] = inftau_PFC_Pyra_Proximal_mh_Na(y(2),Dopamine_Ratio);
y(8) = infh_PFC_Pyra_Proximal_Na;

[tauq_PFC_Pyra_Proximal_KS, infq_PFC_Pyra_Proximal_KS,...
    taur_PFC_Pyra_Proximal_KS, infr_PFC_Pyra_Proximal_KS] = inftau_PFC_Pyra_Proximal_qr_KS(y(2));

y(9) = infq_PFC_Pyra_Proximal_KS;
y(10) = infr_PFC_Pyra_Proximal_KS;

[taua_PFC_Pyra_Distal_A, infa_PFC_Pyra_Distal_A,...
    taub_PFC_Pyra_Distal_A, infb_PFC_Pyra_Distal_A] = inftau_PFC_Pyra_Distal_ab_A(y(3));

y(11) = infa_PFC_Pyra_Distal_A;
y(12) = infb_PFC_Pyra_Distal_A;
y0 = y;




end