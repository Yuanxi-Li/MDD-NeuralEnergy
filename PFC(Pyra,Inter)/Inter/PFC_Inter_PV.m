function [t,y] = PFC_Inter_PV(I_Ext_PFC_Inter_PV)

%%%%%%%%%%%%somatic reduced model
tic
t = 0:0.02:2000;
% I_Ext_Soma = [2:0.2:4.0]
Na_PFC_Inter_PV = PFC_Inter_ionmodel('Na', 35, 5, 55);
K_PFC_Inter_PV = PFC_Inter_ionmodel('K', 9, 5, -90);
Leak_PFC_Inter_PV = PFC_Inter_ionmodel('Leak', 0.1, 1, -65);
V_PFC_Inter_PV = -64;
y0 = Initial_PFC_Inter_PV(V_PFC_Inter_PV);
[t,y] = ode23(@S_PFC_Inter_PV,t,y0);
plot(t,y(:,1))

toc

    function dy = S_PFC_Inter_PV(t,y)
    dy = zeros(3,1);

    %1:Soma 2:Na,h 3:K,n

    %----------------Na 2 h--------------------
    [infm_PFC_Inter_PV_Na,...
        tauh_PFC_Inter_PV_Na, infh_PFC_Inter_PV_Na] = inftau_PFC_Inter_PV_mh_Na(y(1));
    dy(2) = Na_PFC_Inter_PV.phi * (infh_PFC_Inter_PV_Na - y(2)) / tauh_PFC_Inter_PV_Na;
    I_PFC_Inter_PV_Na = Na_PFC_Inter_PV.gmax * (infm_PFC_Inter_PV_Na ^ 3) * y(2) * (y(1) - Na_PFC_Inter_PV.rev);

    %----------------K 2 n--------------------


    [taun_PFC_Inter_PV_K, infn_PFC_Inter_PV_K] = inftau_PFC_Inter_PV_n_K(y(1));
    dy(3) = K_PFC_Inter_PV.phi * (infn_PFC_Inter_PV_K - y(3)) / taun_PFC_Inter_PV_K;
    I_PFC_Inter_PV_K = K_PFC_Inter_PV.gmax * (y(3) ^ 4) * (y(1) - K_PFC_Inter_PV.rev);

    %-----------------Cl--------------------

    I_PFC_Inter_PV_Leak = Leak_PFC_Inter_PV.gmax * (y(1) - Leak_PFC_Inter_PV.rev);

    %Total
    I_PFC_Inter_PV_Ion = I_PFC_Inter_PV_Na+I_PFC_Inter_PV_K+I_PFC_Inter_PV_Leak;
    dy(1) = - I_PFC_Inter_PV_Ion + I_Ext_PFC_Inter_PV;
    end
end