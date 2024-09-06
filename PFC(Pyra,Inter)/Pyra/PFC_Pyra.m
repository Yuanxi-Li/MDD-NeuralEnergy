function [t,y] = PFC_Pyra(I_Ext_PFC_Pyra_Soma,I_Ext_PFC_Pyra_Proximal,I_Ext_PFC_Pyra_Distal)

%%%%%%%%%%%%somatic reduced model
tic
t = 0:0.02:1000;
% I_Ext_Soma = [2:0.2:4.0]
V_PFC_Pyra_Soma = -64.8;
V_PFC_Pyra_Proximal = -64;
V_PFC_Pyra_Distal = -64;
% Name, gmax, phi, rev
Na_PFC_Pyra_Soma = PFC_Pyra_ionmodel('Na', 55, 4, 55);
K_PFC_Pyra_Soma = PFC_Pyra_ionmodel('K', 60, 4, -80);
Ca_PFC_Pyra_Soma = PFC_Pyra_ionmodel('Ca', 1.5, 1, 120);
CaN_PFC_Pyra_Soma = PFC_Pyra_ionmodel('CaN', 0.025, 1, -20);
Leak_PFC_Pyra = PFC_Pyra_ionmodel('Leak', 0.05, 1, -70);
NaP_PFC_Pyra_Proximal = PFC_Pyra_ionmodel('NaP', 0.15, 1, 55);
KS_PFC_Pyra_Proximal = PFC_Pyra_ionmodel('KS', 16, 1, -80); % g = 16 burst
A_PFC_Pyra_Distal = PFC_Pyra_ionmodel('A', 1, 1, -80);
Ca_PFC_Pyra_Distal = PFC_Pyra_ionmodel('Ca', 0.25, 1, 120);


y0 = Initial_PFC_Pyra(V_PFC_Pyra_Soma,V_PFC_Pyra_Proximal,V_PFC_Pyra_Distal);
[t,y] = ode23(@S_PFC_Pyra,t,y0);

toc

    function dy = S_PFC_Pyra(t,y)
        dy = zeros(12,1);  %1 Soma 2 Proximal 3 Distal
                         %4-7 soma 8-10 Proximal 11-12 Distal
        
        %-----------Soma----------------
        %-----------------Soma:Na 4 h--------------------
        [infm_PFC_Pyra_Soma_Na,...
            tauh_PFC_Pyra_Soma_Na, infh_PFC_Pyra_Soma_Na] = inftau_PFC_Pyra_Soma_mh_Na(y(1));

        dy(4) = Na_PFC_Pyra_Soma.phi*(infh_PFC_Pyra_Soma_Na-y(4))/tauh_PFC_Pyra_Soma_Na;
        I_PFC_Pyra_Soma_Na = Na_PFC_Pyra_Soma.gmax*(infm_PFC_Pyra_Soma_Na^3)*y(4)*(y(1)-Na_PFC_Pyra_Soma.rev);
        
        %-----------------Soma:K 5 n--------------------
        
        [taun_PFC_Pyra_Soma_K, infn_PFC_Pyra_Soma_K] = inftau_PFC_Pyra_Soma_n_K(y(1));
        dy(5) = K_PFC_Pyra_Soma.phi*(infn_PFC_Pyra_Soma_K-y(5))/taun_PFC_Pyra_Soma_K;
        I_PFC_Pyra_Soma_K = K_PFC_Pyra_Soma.gmax*(y(5)^4)*(y(1)-K_PFC_Pyra_Soma.rev);
        
        %-----------------Soma:Ca m--------------------
        
        infm_PFC_Pyra_Soma_Ca = 1/(1+exp(-(y(1)+20)/9));
        I_PFC_Pyra_Soma_Ca = Ca_PFC_Pyra_Soma.gmax*(infm_PFC_Pyra_Soma_Ca^2)*(y(1)-Ca_PFC_Pyra_Soma.rev);
        
        %-----------------Soma:CaN m 6--------------------
        
        alpham_PFC_Pyra_Soma_CaN = 0.0056;
        betam_PFC_Pyra_Soma_CaN = 0.002;
        infm_PFC_Pyra_Soma_CaN = alpham_PFC_Pyra_Soma_CaN*(y(7)^2)/(alpham_PFC_Pyra_Soma_CaN*(y(7)^2)+betam_PFC_Pyra_Soma_CaN);
        taum_PFC_Pyra_Soma_CaN = 1/(alpham_PFC_Pyra_Soma_CaN*(y(7)^2)+betam_PFC_Pyra_Soma_CaN);
        dy(6) = (infm_PFC_Pyra_Soma_CaN-y(6))/taum_PFC_Pyra_Soma_CaN;
        I_PFC_Pyra_Soma_CaN = CaN_PFC_Pyra_Soma.gmax*(y(6)^2)*(y(1)-CaN_PFC_Pyra_Soma.rev);
        
        %-----------------Soma:Ca dynamics 7--------------------
        
        alpha_PFC_Pyra_Soma_Cadyn = 0.000667;
        tau_PFC_Pyra_Soma_Cadyn = 240;
        I_PFC_Pyra_Soma_Cadyn = I_PFC_Pyra_Soma_Ca+I_PFC_Pyra_Soma_CaN;
        dy(7) = -alpha_PFC_Pyra_Soma_Cadyn*I_PFC_Pyra_Soma_Cadyn-y(7)/tau_PFC_Pyra_Soma_Cadyn;
        
        %-----------------Soma:Cl--------------------
        I_PFC_Pyra_Soma_Leak = Leak_PFC_Pyra.gmax*(y(1)-Leak_PFC_Pyra.rev);
        
        
        
        %-----------Proximal Dendrites----------------
        %-----------------Proximal:NaP 8 h--------------------
        
        [infm_PFC_Pyra_Proximal_Na,...
    tauh_PFC_Pyra_Proximal_Na, infh_PFC_Pyra_Proximal_Na] = inftau_PFC_Pyra_Proximal_mh_Na(y(2));
        dy(8) = (infh_PFC_Pyra_Proximal_Na-y(8))/tauh_PFC_Pyra_Proximal_Na;
        I_PFC_Pyra_Proximal_NaP = NaP_PFC_Pyra_Proximal.gmax*(infm_PFC_Pyra_Proximal_Na^3)*y(8)*(y(2)-NaP_PFC_Pyra_Proximal.rev);
        
        %-----------------Proximal:KS q r 9 10--------------------
        
        %        16;
        [tauq_PFC_Pyra_Proximal_KS, infq_PFC_Pyra_Proximal_KS,...
            taur_PFC_Pyra_Proximal_KS, infr_PFC_Pyra_Proximal_KS] = inftau_PFC_Pyra_Proximal_qr_KS(y(2));
        dy(9) = (infq_PFC_Pyra_Proximal_KS-y(9))/tauq_PFC_Pyra_Proximal_KS;
        dy(10) = (infr_PFC_Pyra_Proximal_KS-y(10))/taur_PFC_Pyra_Proximal_KS;
        I_PFC_Pyra_Proximal_KS = KS_PFC_Pyra_Proximal.gmax*y(9)*y(10)*(y(2)-KS_PFC_Pyra_Proximal.rev);
        
        
        %-----------------Proximal:Cl--------------------

        I_PFC_Pyra_Proximal_Leak =  Leak_PFC_Pyra.gmax*(y(2)-Leak_PFC_Pyra.rev);
        
        
        %-----------Distal Dendrites----------------
        %-----------------Distal:A a b 11 12--------------------
        
        
        [taua_PFC_Pyra_Distal_A, infa_PFC_Pyra_Distal_A,...
            taub_PFC_Pyra_Distal_A, infb_PFC_Pyra_Distal_A] = inftau_PFC_Pyra_Distal_ab_A(y(3));
        dy(11) = (infa_PFC_Pyra_Distal_A-y(11))/taua_PFC_Pyra_Distal_A;
        dy(12) = (infb_PFC_Pyra_Distal_A-y(12))/taub_PFC_Pyra_Distal_A;
        I_PFC_Pyra_Distal_A = A_PFC_Pyra_Distal.gmax*(y(11)^4)*y(12)*(y(3)-A_PFC_Pyra_Distal.rev);
        
        %-----------------Distal:Ca--------------------
        infm_PFC_Pyra_Distal_Ca = 1/(1+exp(-(y(3)+20)/9));
        I_PFC_Pyra_Distal_Ca = Ca_PFC_Pyra_Distal.gmax*(infm_PFC_Pyra_Distal_Ca^2)*(y(3)-Ca_PFC_Pyra_Distal.rev);        
        
        
        %-----------------Distal:Cl--------------------

        I_PFC_Pyra_Distal_Leak =  Leak_PFC_Pyra.gmax*(y(3)-Leak_PFC_Pyra.rev);

        
        
        %-----------------Total--------------------
       
        I_PFC_Pyra_Soma_Ion = I_PFC_Pyra_Soma_Na+I_PFC_Pyra_Soma_K+I_PFC_Pyra_Soma_Ca+I_PFC_Pyra_Soma_CaN+I_PFC_Pyra_Soma_Leak;
%         I_Soma = I_Na_Soma+I_K_Soma;
        gc1 = 0.75;
        p1 = 0.5;
        dy(1) = -I_PFC_Pyra_Soma_Ion-gc1*(y(1)-y(2))/p1+I_Ext_PFC_Pyra_Soma;
        I_PFC_Pyra_Proximal_Ion = I_PFC_Pyra_Proximal_NaP+I_PFC_Pyra_Proximal_KS+I_PFC_Pyra_Proximal_Leak;
        gc2 = 0.25;
        p2 = 0.3;
        dy(2) = -I_PFC_Pyra_Proximal_Ion-gc1*(y(2)-y(1))/p2 -gc2*(y(2)-y(3))/p2+I_Ext_PFC_Pyra_Proximal;
        I_PFC_Pyra_Distal_Ion = I_PFC_Pyra_Distal_A+I_PFC_Pyra_Distal_Ca+I_PFC_Pyra_Distal_Leak;
%         I_Distal = I_A_Distal+I_L_Distal;
        dy(3) = -I_PFC_Pyra_Distal_Ion-gc2*(y(3)-y(2))/(1-p1-p2)+I_Ext_PFC_Pyra_Distal; 
        
    end
end