function [infm_PFC_Pyra_Proximal_Na,...
    tauh_PFC_Pyra_Proximal_Na, infh_PFC_Pyra_Proximal_Na] = inftau_PFC_Pyra_Proximal_mh_Na(Vm,Dopamine_Ratio)
% infm_PFC_Pyra_Proximal_Na = 1/(1+exp(-(Vm+55.7)/7.7)); % Standard
% alphah_PFC_Pyra_Proximal_Na = 0.001*exp((-85-Vm)/30); % Standard
% betah_PFC_Pyra_Proximal_Na = 0.0034/(exp((-17-Vm)/10)+1);% Standard
% infm_PFC_Pyra_Proximal_Na = 1/(1+exp(-(Vm+60.7)/7.7));%Dopamine
% alphah_PFC_Pyra_Proximal_Na = 0.0005*exp((-85-Vm)/30);%Dopamine
% betah_PFC_Pyra_Proximal_Na = 0.0017/(exp((-17-Vm)/10)+1);%Dopamine
infm_PFC_Pyra_Proximal_Na = 1/(1+exp(-(Vm+(55.7+5*Dopamine_Ratio))/7.7));%Dopamine Ratio
alphah_PFC_Pyra_Proximal_Na = (0.001-0.0005*Dopamine_Ratio)*exp((-85-Vm)/30);%Dopamine Ratio
betah_PFC_Pyra_Proximal_Na = (0.0034-0.0017*Dopamine_Ratio)/(exp((-17-Vm)/10)+1);%Dopamine Ratio
[infh_PFC_Pyra_Proximal_Na,tauh_PFC_Pyra_Proximal_Na] = alpha_beta_solve(alphah_PFC_Pyra_Proximal_Na, betah_PFC_Pyra_Proximal_Na);
end