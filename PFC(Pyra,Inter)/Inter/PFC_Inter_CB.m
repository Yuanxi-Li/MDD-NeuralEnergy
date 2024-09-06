function [t,y] = PFC_Inter_CB(I_Ext_PFC_Inter_CB)

%%%%%%%%%%%%somatic reduced model
tic
t = 0:0.02:2000;
% I_Ext_Soma = [2:0.2:4.0]
V_PFC_Inter_CB = -64;
y0 = Initial_PFC_Inter_CB(V_PFC_Inter_CB);
Na_PFC_Inter_CB = PFC_Inter_ionmodel('Na', 35, 5, 55);
K_PFC_Inter_CB = PFC_Inter_ionmodel('K', 9, 5, -85);
Ca_PFC_Inter_CB = PFC_Pyra_ionmodel('Ca', 1, 1, 120);
KCa_PFC_Inter_CB = PFC_Pyra_ionmodel('KCa', 1, 1, -85);
h_PFC_Inter_CB = PFC_Pyra_ionmodel('KCa', 0.15, 1, -40);
Leak_PFC_Inter_CB = PFC_Pyra_ionmodel('Leak', 0.1, 1, -65);


[t,y] = ode23(@S_PFC_Inter_CB,t,y0);

toc

    function dy = S_PFC_Inter_CB(t,y)
    dy = zeros(5,1);


    %'''----------------Na 2 h--------------------'''
    
    [infm_PFC_Inter_CB_Na,...
        tauh_PFC_Inter_CB_Na, infh_PFC_Inter_CB_Na] = inftau_PFC_Inter_CB_mh_Na(y(1));
    dy(2) = Na_PFC_Inter_CB.phi * (infh_PFC_Inter_CB_Na - y(2)) / tauh_PFC_Inter_CB_Na;
    I_PFC_Inter_CB_Na = Na_PFC_Inter_CB.gmax * (infm_PFC_Inter_CB_Na ^ 3) * y(2) * (y(1) - Na_PFC_Inter_CB.rev);

    %'''----------------K 3 n--------------------'''

    [taun_PFC_Inter_CB_K, infn_PFC_Inter_CB_K] = inftau_PFC_Inter_CB_n_K(y(1));
    dy(3) = K_PFC_Inter_CB.phi * (infn_PFC_Inter_CB_K - y(3)) / taun_PFC_Inter_CB_K;
    I_PFC_Inter_CB_K = K_PFC_Inter_CB.gmax * (y(3) ^ 4) * (y(1) - K_PFC_Inter_CB.rev);

    %'''----------------Ca m--------------------'''
    infm_PFC_Inter_CB_Ca = 1 / (1 + exp(-(y(1) + 20) / 9));
    I_PFC_Inter_CB_Ca = Ca_PFC_Inter_CB.gmax * (infm_PFC_Inter_CB_Ca ^ 2) * (y(1) - Ca_PFC_Inter_CB.rev);


    %'''-----------------Ca dynamics 4--------------------'''
    alpha_Cadyn_Inter_CB = 0.002;
    tau_Cadyn_Inter_CB = 80;
    I_PFC_Inter_CB_Cadyn = I_PFC_Inter_CB_Ca;
    dy(4) = -alpha_Cadyn_Inter_CB * I_PFC_Inter_CB_Cadyn - y(4) / tau_Cadyn_Inter_CB;

    %'''-----------------KCa--------------------'''
    KD_Inter_CB = 30;
    I_PFC_Inter_CB_KCa = KCa_PFC_Inter_CB.gmax * (y(4) / (y(4) + KD_Inter_CB)) * (y(1) - KCa_PFC_Inter_CB.rev);

    %'''-----------------Ih H 5--------------------'''

    [tauh_PFC_Inter_CB_h, infh_PFC_Inter_CB_h] = inftau_PFC_Inter_CB_h_h(y(1));
    dy(5) = (infh_PFC_Inter_CB_h - y(5)) / tauh_PFC_Inter_CB_h;
    I_PFC_Inter_CB_h = h_PFC_Inter_CB.gmax * y(5) * (y(1) - h_PFC_Inter_CB.rev);

    %'''-----------------Cl--------------------'''

    I_PFC_Inter_CB_Leak = Leak_PFC_Inter_CB.gmax * (y(1) - Leak_PFC_Inter_CB.rev);

    %'''Total'''
    I_PFC_Inter_CB_Ion = I_PFC_Inter_CB_Na + I_PFC_Inter_CB_K + I_PFC_Inter_CB_Ca + ... 
        I_PFC_Inter_CB_KCa + I_PFC_Inter_CB_h + I_PFC_Inter_CB_Leak;
    dy(1) = -I_PFC_Inter_CB_Ion + I_Ext_PFC_Inter_CB;
    
    end
end