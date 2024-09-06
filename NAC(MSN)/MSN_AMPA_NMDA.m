function [t,y] = MSN_AMPA_NMDA(I_Ext_MSN_Soma,I_Ext_MSN_Dendrite)
%%%%%%%%% S = 2.7*1e-5 (cm^2);
tic
%%%%%%%depth = 15;
[t1,y1] = PFC_Pyra(3,3,3);
ypre = y1(:,1);
t = 0:0.02:1001;
% I_Ext = I_Clamp+2*sin(0.05*t);
V_MSN_Soma = -90;
V_MSN_Dendrite = -90;
y0 = Initial_MSN(V_MSN_Soma,V_MSN_Dendrite);
y0 = [y0;0;0];
% Parameters
NaF_MSN_Soma = MSN_ionmodel1('NaF', 1.5, nan, -23.9, -62.9, -11.8, 10.7, 50);
NaF_MSN_Dendrite = MSN_ionmodel1('NaF', 0.0195, nan, -23.9, -62.9, -11.8, 10.7, 50);
NaP_MSN_Soma = MSN_ionmodel1('NaP', 4e-5, nan, -52.6, -48.8, -4.6, 10, 50);
NaP_MSN_Dendrite = MSN_ionmodel1('NaP', 1.38e-7, nan, -52.6, -48.8, -4.6, 10, 50);
KAf_MSN_Soma = MSN_ionmodel1('KAf', 0.225, nan, -10, -75.6, -17.7, 10, -90);
KAf_MSN_Dendrite = MSN_ionmodel1('KAf', 0.225, nan, -10, -75.6, -17.7, 10, -90);
KAs_MSN_Soma = MSN_ionmodel1('KAs', 0.0104, 0.996, -27, -33.5, -16, 21.5, -90);
KAs_MSN_Dendrite = MSN_ionmodel1('KAs', 0.0104, 0.996, -27, -33.5, -16, 21.5, -90);
KIR_MSN = MSN_ionmodel1('KIR', 1.4e-4, nan, -82, nan, 13, nan, -90);
KRP_MSN = MSN_ionmodel1('KRP', 0.001, 0.7, -13.5, -54.7, -11.8, 18.6, -90);
Leak_MSN = MSN_ionmodel1('Leak', 11.5e-6, nan, nan, nan, nan, nan, -70);
BKKCa_MSN = MSN_ionmodel1('BKKCa', 0.001, nan, nan, nan, nan, nan, -90);
SKKCa_MSN = MSN_ionmodel1('SKKCa', 0.145, nan, nan, nan, nan, nan, -90);

CaL1_2_MSN = MSN_ionmodel2('CaL1.2', 6.7e-6, 0.17, -8.9, -13.4, -6.7, 11.9);
CaL1_3_MSN = MSN_ionmodel2('CaL1.3', 4.25e-7, nan, -33, -13.4, -6.7, 11.9);
CaN_MSN = MSN_ionmodel2('CaN', 1.0e-5, 0.21, -8.7, -74.8, -7.4, 6.5);
CaQ_MSN = MSN_ionmodel2('CaQ', 6.0e-6, nan, -9.0, nan, -6.6, nan);
CaR_MSN = MSN_ionmodel2('CaR', 2.6e-5, nan, -10.3, -33.3, -6.6, 17);
CaT_MSN = MSN_ionmodel2('CaT', 4e-7, nan, -51.73, -80, -6.53, 6.7);

g_Rall_MSN = MSN_RallModel;
AMPA_PFCPyra2MSN = SynapseChannel('AMPA','PFCPyra','MSN',1,1,0.2,0);
NMDA_PFCPyra2MSN = SynapseChannel('NMDA','PFCPyra','MSN',1,0.072,0.0067,0);


[t,y] = ode45(@S_MSN_AMPA_NMDA,t,y0);
plot(t,y(:,1))
toc


    function dy = S_MSN_AMPA_NMDA(t,y)
        dy = zeros(58,1);%---------------Column      1 Soma ; 2 Proximal
                     %---------------           Soma: 3 - 30
                     %---------------           3,4 NaF m h ; 5,6 NaP m h ; 7,8 KAf m h;
                     %---------------           9,10 KAs m h ; 11 KIR m ; 12,13 KRP m h ; 14,15 CaL1_2 m h ;  
                     %---------------           16,17 CaL1_3 m h ; 18,19 CaN m h ; 20 CaQ m ; 21,22 CaR m h ;  
                     %---------------           23,24 CaT m h; 
                     %---------------           25 L型Ca CaL_i  ; 26:Ca_i
                     %---------------           27 SKKCa
                     %---------------           28,29,30:BKKCa    28:open,29:close,30:inactive
                     %---------------           Proximal: 31-56,no KRP
                     %---------------           31,32 NaF m h ; 33,34 NaP m h ; 35,36 KAf m h;
                     %---------------           37,38 KAs m h ; 39 KIR m ; 40,41 CaL1_2 m h ;  
                     %---------------           42,43 CaL1_3 m h ; 44,45 CaN m h ; 46 CaQ m ; 47,48 CaR m h ;  
                     %---------------           49,50 CaT m h; 
                     %---------------           51 L型Ca CaL_i  ; 52:Ca_i
                     %---------------           53 SKKCa
                     %---------------           54,55,56:BKKCa    54:open,55:close,56:inactive
        
        %----------------------NaF----------------------------
        
        %计算NaF的mtau htau
        [tauh_NaF_MSN_Soma, taum_NaF_MSN_Soma] = tau_NaF_MSN(y(1));

        %计算mInf hInf
        mInf_NaF_MSN_Soma = infinite(y(1),NaF_MSN_Soma.mVhalf,NaF_MSN_Soma.mk);
        hInf_NaF_MSN_Soma = infinite(y(1),NaF_MSN_Soma.hVhalf,NaF_MSN_Soma.hk);
        
        %计算3 4m,h
        dy(3) = (mInf_NaF_MSN_Soma-y(3))/taum_NaF_MSN_Soma;
        dy(4) = (hInf_NaF_MSN_Soma-y(4))/tauh_NaF_MSN_Soma;
        I_NaF_MSN_Soma = NaF_MSN_Soma.gmax*(y(3)^3)*y(4)*(y(1)-NaF_MSN_Soma.rev);%-----mA
        
        %----------------------NaP----------------------------
        
        %计算NaP的mtau htau
        [tauh_NaP_MSN_Soma, taum_NaP_MSN_Soma] = tau_NaP_MSN(y(1));
        
        %计算mInf hInf
        mInf_NaP_MSN_Soma = infinite(y(1),NaP_MSN_Soma.mVhalf,NaP_MSN_Soma.mk);
        hInf_NaP_MSN_Soma = infinite(y(1),NaP_MSN_Soma.hVhalf,NaP_MSN_Soma.hk);
        
        %计算5 6  m,h
        dy(5) = (mInf_NaP_MSN_Soma-y(5))/taum_NaP_MSN_Soma;
        dy(6) = (hInf_NaP_MSN_Soma-y(6))/tauh_NaP_MSN_Soma;
        I_NaP_MSN_Soma = NaP_MSN_Soma.gmax*y(5)*y(6)*(y(1)-NaP_MSN_Soma.rev);%---------mA
        
        %----------------------KAf----------------------------
        
        %计算KAf的mtau htau
        [tauh_KAf_MSN_Soma, taum_KAf_MSN_Soma] = tau_KAf_MSN(y(1));
        
        %计算mInf hInf
        mInf_KAf_MSN_Soma = infinite(y(1),KAf_MSN_Soma.mVhalf,KAf_MSN_Soma.mk);
        hInf_KAf_MSN_Soma = infinite(y(1),KAf_MSN_Soma.hVhalf,KAf_MSN_Soma.hk);
        
        %计算7 8m,h
        dy(7) = (mInf_KAf_MSN_Soma-y(7))/taum_KAf_MSN_Soma;
        dy(8) = (hInf_KAf_MSN_Soma-y(8))/tauh_KAf_MSN_Soma;
        I_KAf_MSN_Soma = KAf_MSN_Soma.gmax*(y(7)^2)*y(8)*(y(1)-KAf_MSN_Soma.rev);%---------mA
                
        
        %----------------------KAs----------------------------
       
        
        %计算KAs的mtau htau
        [tauh_KAs_MSN_Soma, taum_KAs_MSN_Soma] = tau_KAs_MSN(y(1));
        
        %计算mInf hInf
        mInf_KAs_MSN_Soma = infinite(y(1),KAs_MSN_Soma.mVhalf,KAs_MSN_Soma.mk);
        hInf_KAs_MSN_Soma = infinite(y(1),KAs_MSN_Soma.hVhalf,KAs_MSN_Soma.hk);
        
        %计算9 10m,h
        dy(9) = (mInf_KAs_MSN_Soma-y(9))/taum_KAs_MSN_Soma;
        dy(10) = (hInf_KAs_MSN_Soma-y(10))/tauh_KAs_MSN_Soma;
        I_KAs_MSN_Soma = KAs_MSN_Soma.gmax*(y(9)^2)*(KAs_MSN_Soma.a*y(10)+(1-KAs_MSN_Soma.a))*(y(1)-KAs_MSN_Soma.rev);%---------mA
        
        %----------------------KIR----------------------------

        
        %计算KIR的mtau
        taum_KIR_MSN_Soma = tau_KIR_MSN(y(1));
        
        %计算mInf
        mInf_KIR_MSN_Soma = infinite(y(1),KIR_MSN.mVhalf,KIR_MSN.mk);
        
        %计算11 m
        dy(11) = (mInf_KIR_MSN_Soma-y(11))/taum_KIR_MSN_Soma;
        I_KIR_MSN_Soma = KIR_MSN.gmax*y(11)*(y(1)-KIR_MSN.rev);%---------mA
        
        %----------------------KRP----------------------------
        
        %计算KRP的mtau htau
        [tauh_KRP_MSN_Soma, taum_KRP_MSN_Soma] = tau_KRP_MSN(y(1));
        
        %计算mInf hInf
        mInf_KRP_MSN_Soma = infinite(y(1),KRP_MSN.mVhalf,KRP_MSN.mk);
        hInf_KRP_MSN_Soma = infinite(y(1),KRP_MSN.mVhalf,KRP_MSN.mk);
        
        %计算12 13 m,h
        dy(12) = (mInf_KRP_MSN_Soma-y(12))/taum_KRP_MSN_Soma;
        dy(13) = (hInf_KRP_MSN_Soma-y(13))/tauh_KRP_MSN_Soma;
        I_KRP_MSN_Soma = KRP_MSN.gmax*y(12)*(KRP_MSN.a*y(13)+(1-KRP_MSN.a))*(y(1)-KRP_MSN.rev);%---------mA
        
        %----------------------Ca固定参数----------------------------
        
        z_MSN = 2;
        F_MSN = 96489;%-------------Faraday(C/mod)
        R_MSN = 8.31;%------------J/(mol*K)
        T_MSN = 35+273.15;%--------K
        Ca_0_MSN = 5;%---------Ca库
        CaL_0_MSN = 5;%--------CaL库
        
        

        
        %----------------------CaL1_2----------------------------
        
        %计算CaL1_2的mtau htau
        [tauh_CaL_MSN_Soma, taum_CaL_MSN_Soma] = tau_CaL_MSN(y(1));
        
        %计算mInf hInf
        mInf_CaL1_2_MSN_Soma = infinite(y(1),CaL1_2_MSN.mVhalf,CaL1_2_MSN.mk);
        hInf_CaL1_2_MSN_Soma = infinite(y(1),CaL1_2_MSN.mVhalf,CaL1_2_MSN.mk);
        
        %计算 14 15 m,h  % 25 L型Ca CaL_i
        dy(14) = (mInf_CaL1_2_MSN_Soma-y(14))/taum_CaL_MSN_Soma;
        dy(15) = (hInf_CaL1_2_MSN_Soma-y(15))/tauh_CaL_MSN_Soma;
        
        g_CaL1_2_MSN_Soma = g_Ca_MSN(y(1),y(25),CaL_0_MSN);
        I_CaL1_2_MSN_Soma = g_CaL1_2_MSN_Soma*CaL1_2_MSN.Pbar*y(14)*y(14)*(y(15)*CaL1_2_MSN.a+(1-CaL1_2_MSN.a));%------mA
        
        %----------------------CaL1_3----------------------------
        
        %计算CaL1_3的mtau htau-----------与CaL1_2的相同
        
        %计算mInf hInf
        mInf_CaL1_3_MSN_Soma = infinite(y(1),CaL1_3_MSN.mVhalf,CaL1_3_MSN.mk);
        hInf_CaL1_3_MSN_Soma = infinite(y(1),CaL1_3_MSN.mVhalf,CaL1_3_MSN.mk);
        
        %计算16 17m,h  %25 L型Ca CaL_i
        dy(16) = (mInf_CaL1_3_MSN_Soma-y(16))/taum_CaL_MSN_Soma;
        dy(17) = (hInf_CaL1_3_MSN_Soma-y(17))/tauh_CaL_MSN_Soma;
        
        g_CaL1_3_MSN_Soma = g_Ca_MSN(y(1),y(25),CaL_0_MSN);
        I_CaL1_3_MSN_Soma = g_CaL1_3_MSN_Soma*CaL1_3_MSN.Pbar*y(16)*y(16)*y(17);%------mA

        
         
        %----------------------CaN----------------------------
        
        
        %计算CaN的mtau htau
        [tauh_CaN_MSN_Soma, taum_CaN_MSN_Soma] = tau_CaN_MSN(y(1));


        
        %计算mInf hInf
        mInf_CaN_MSN_Soma = infinite(y(1), CaN_MSN.mVhalf, CaN_MSN.mk);
        hInf_CaN_MSN_Soma = infinite(y(1), CaN_MSN.hVhalf, CaN_MSN.hk);
        
        %计算 18 19 m,h  %26 Ca Ca_i
        dy(18) = (mInf_CaN_MSN_Soma-y(18))/taum_CaN_MSN_Soma;
        dy(19) = (hInf_CaN_MSN_Soma-y(19))/tauh_CaN_MSN_Soma;
        
        g_CaN_MSN_Soma = g_Ca_MSN(y(1),y(26),Ca_0_MSN);
        I_CaN_MSN_Soma = g_CaN_MSN_Soma*CaN_MSN.Pbar*y(18)*y(18)*(CaN_MSN.a*y(19)+(1-CaN_MSN.a));%------mA

             
        %----------------------CaQ----------------------------
        
        
        %计算CaQ的mtau
        taum_CaQ_MSN_Soma = tau_CaQ_MSN(y(1));
        
        %计算mInf hInf
        mInf_CaQ_MSN_Soma = infinite(y(1),CaQ_MSN.mVhalf,CaQ_MSN.mk);

        
        %计算 20 m  % 26 Ca Ca_i
        dy(20) = (mInf_CaQ_MSN_Soma-y(20))/taum_CaQ_MSN_Soma;
        
        g_CaQ_MSN_Soma = g_Ca_MSN(y(1),y(26),Ca_0_MSN);
        I_CaQ_MSN_Soma = g_CaQ_MSN_Soma*CaQ_MSN.Pbar*y(20)*y(20);%------mA
        
        
        
        %----------------------CaR----------------------------
        

        
        %计算CaR的mtau htau
        [tauh_CaR_MSN_Soma, taum_CaR_MSN_Soma] = tau_CaR_MSN(y(1));
        
        %计算mInf hInf
        mInf_CaR_MSN_Soma = infinite(y(1),CaR_MSN.mVhalf,CaR_MSN.mk);
        hInf_CaR_MSN_Soma = infinite(y(1),CaR_MSN.hVhalf,CaR_MSN.hk);
        
        %计算 21 22 m,h  % 26 Ca Ca_i
        dy(21) = (mInf_CaR_MSN_Soma-y(21))/taum_CaR_MSN_Soma;
        dy(22) = (hInf_CaR_MSN_Soma-y(22))/tauh_CaR_MSN_Soma;
        
        g_CaR_MSN_Soma = g_Ca_MSN(y(1),y(26),Ca_0_MSN);
        I_CaR_MSN_Soma = g_CaR_MSN_Soma*CaR_MSN.Pbar*y(21)*y(21)*y(21)*y(22);%------mA
         
        
        %----------------------CaT----------------------------
        
        %计算CaT的mtau htau
        [tauh_CaT_MSN_Soma, taum_CaT_MSN_Soma] = tau_CaT_MSN(y(1));
        %计算mInf hInf
        mInf_CaT_MSN_Soma = infinite(y(1),CaT_MSN.mVhalf,CaT_MSN.mk);
        hInf_CaT_MSN_Soma = infinite(y(1),CaT_MSN.hVhalf,CaT_MSN.hk);
        
        %计算 23 24 m,h  % 25 CaL CaL_i
        dy(23) = (mInf_CaT_MSN_Soma-y(23))/taum_CaT_MSN_Soma;
        dy(24) = (hInf_CaT_MSN_Soma-y(24))/tauh_CaT_MSN_Soma;
        
        g_CaT_MSN_Soma = g_Ca_MSN(y(1),y(25),Ca_0_MSN);
        I_CaT_MSN_Soma = g_CaT_MSN_Soma*CaT_MSN.Pbar*y(23)*y(23)*y(23)*y(24);%------mA
 
        %----------------------25 CaLdyn ----------------------------        

        depth_MSN = 0.1;%--------Shell depth,miumeter
        drive_MSN = 10000;%-------influx
        pump_MSN = 0.02;%-------efflux
        Kt_Ca_MSN = 1e-4;%-------mmol/L/ms
        Ca_i_inf_MSN = 1e-5;%------mmol/L
        CaL_i_inf_MSN = 1e-5;
        Kd_Ca_MSN = 1e-4;%--------mmol/L
        tau_R_MSN = 43;
        I_CaL_MSN_Soma = I_CaL1_2_MSN_Soma+I_CaL1_3_MSN_Soma+I_CaT_MSN_Soma;
        drive_channel_CaL_MSN_Soma = -drive_MSN * I_CaL_MSN_Soma /(2*F_MSN*depth_MSN);%----------把传入的钙（从通道）转换成内部的浓度变化
        %---------------------如果<=0,不能泵入
        if drive_channel_CaL_MSN_Soma<=0
            drive_channel_CaL_MSN_Soma = 0;
        end
        drive_pump_CaL_MSN_Soma = -Kt_Ca_MSN * y(25)/(y(25)+Kd_Ca_MSN);
        dy(25) = drive_channel_CaL_MSN_Soma + pump_MSN*drive_pump_CaL_MSN_Soma +(CaL_i_inf_MSN-y(25))/tau_R_MSN;
        
        
        %----------------------26 Ca ----------------------------        
        I_Ca_MSN_Soma = I_CaN_MSN_Soma+I_CaQ_MSN_Soma+I_CaR_MSN_Soma;
        drive_channel_Ca_MSN_Soma = -drive_MSN * I_Ca_MSN_Soma /(2*F_MSN*depth_MSN);%----------把传入的钙（从通道）转换成内部的浓度变化
        if drive_channel_Ca_MSN_Soma<=0
            drive_channel_Ca_MSN_Soma = 0;
        end
        drive_pump_Ca_MSN_Soma = -Kt_Ca_MSN * y(26)/(y(26)+Kd_Ca_MSN);
        dy(26) = drive_channel_Ca_MSN_Soma + pump_MSN*drive_pump_Ca_MSN_Soma +(Ca_i_inf_MSN-y(26))/tau_R_MSN;
        
        %----------------------27 SKKCa ----------------------------
        [tau_SKKCa_MSN_Soma,oinf_SKKCa_MSN_Soma] = rate_SKKCa_MSN(y(1),y(26));
        dy(27) = (oinf_SKKCa_MSN_Soma-y(27))/tau_SKKCa_MSN_Soma;
        I_SKKCa_MSN_Soma = SKKCa_MSN.gmax*y(27)*(y(1)-SKKCa_MSN.rev);

        
        %----------------------BKKCa ----------------------------
        %---------------28 29 30 :BKKCa    28:open,29:close,30:inactive
        [k1_BKKCa_MSN_Soma,k2_BKKCa_MSN_Soma,k3_BKKCa_MSN_Soma,k4_BKKCa_MSN_Soma] = rates_BKKCa_MSN(y(1),y(26));
        dy(28) = -k4_BKKCa_MSN_Soma*y(28)-k1_BKKCa_MSN_Soma*y(28)+k3_BKKCa_MSN_Soma*y(29);%%%%%%%%open state
        dy(29) = -k3_BKKCa_MSN_Soma*y(29)+k4_BKKCa_MSN_Soma*y(28)+k2_BKKCa_MSN_Soma*y(30);%%%%%%%%close state
        dy(30) = k1_BKKCa_MSN_Soma*y(28)-k2_BKKCa_MSN_Soma*y(30);%%%%%%%inactive state
        I_BKKCa_MSN_Soma = BKKCa_MSN.gmax*y(28)*(y(1)-BKKCa_MSN.rev);
        
        
        
        %----------------------I_L漏电流----------------------------             
        I_L_MSN_Soma = Leak_MSN.gmax*(y(1)-Leak_MSN.rev);
        
        %----------------------Proximal 2,31-56----------------------------
        %---------------           31,32 NaF m h ; 33,34 NaP m h ; 35,36 KAf m h;
        %---------------           37,38 KAs m h ; 39 KIR m ; 40,41 CaL1_2 m h ;
        %---------------           42,43 CaL1_3 m h ; 44,45 CaN m h ; 46 CaQ m ; 47,48 CaR m h ;
        %---------------           59,50 CaT m h;
        %---------------           51 L型Ca CaL_i  ; 52:Ca_i
        %---------------           53 SKKCa
        %---------------           54,55,56:BKKCa    54:open,55:close,56:inactive
       
        %----------------------NaF---------------------------
        
        %计算NaF的mtau htau
        [tauh_NaF_MSN_Dendrite, taum_NaF_MSN_Dendrite] = tau_NaF_MSN(y(2));

        %计算mInf hInf
        mInf_NaF_MSN_Dendrite = infinite(y(2),NaF_MSN_Dendrite.mVhalf,NaF_MSN_Dendrite.mk);
        hInf_NaF_MSN_Dendrite = infinite(y(2),NaF_MSN_Dendrite.hVhalf,NaF_MSN_Dendrite.hk);
        
        %计算31 32 m,h
        dy(31) = (mInf_NaF_MSN_Dendrite-y(31))/taum_NaF_MSN_Dendrite;
        dy(32) = (hInf_NaF_MSN_Dendrite-y(32))/tauh_NaF_MSN_Dendrite;
        I_NaF_MSN_Dendrite = NaF_MSN_Dendrite.gmax*(y(31)^3)*y(32)*(y(2)-NaF_MSN_Dendrite.rev);%-----mA
        
        %----------------------NaP----------------------------
        
        %计算NaP的mtau htau
        [tauh_NaP_MSN_Dendrite, taum_NaP_MSN_Dendrite] = tau_NaP_MSN(y(2));
        
        %计算mInf hInf
        mInf_NaP_MSN_Dendrite = infinite(y(2),NaP_MSN_Dendrite.mVhalf,NaP_MSN_Dendrite.mk);
        hInf_NaP_MSN_Dendrite = infinite(y(2),NaP_MSN_Dendrite.hVhalf,NaP_MSN_Dendrite.hk);
        
        %计算 33 34  m,h
        dy(33) = (mInf_NaP_MSN_Dendrite-y(33))/taum_NaP_MSN_Dendrite;
        dy(34) = (hInf_NaP_MSN_Dendrite-y(34))/tauh_NaP_MSN_Dendrite;
        I_NaP_MSN_Dendrite = NaP_MSN_Dendrite.gmax*y(33)*y(34)*(y(2)-NaP_MSN_Dendrite.rev);%---------mA
        
        %----------------------KAf----------------------------
        
        %计算KAf的mtau htau
        [tauh_KAf_MSN_Dendrite, taum_KAf_MSN_Dendrite] = tau_KAf_MSN(y(2));
        
        %计算mInf hInf
        mInf_KAf_MSN_Dendrite = infinite(y(2),KAf_MSN_Dendrite.mVhalf,KAf_MSN_Dendrite.mk);
        hInf_KAf_MSN_Dendrite = infinite(y(2),KAf_MSN_Dendrite.hVhalf,KAf_MSN_Dendrite.hk);
        
        %计算35 36 m,h
        dy(35) = (mInf_KAf_MSN_Dendrite-y(35))/taum_KAf_MSN_Dendrite;
        dy(36) = (hInf_KAf_MSN_Dendrite-y(36))/tauh_KAf_MSN_Dendrite;
        I_KAf_MSN_Dendrite = KAf_MSN_Dendrite.gmax*(y(35)^2)*y(36)*(y(2)-KAf_MSN_Dendrite.rev);%---------mA
                
        
        %----------------------KAs----------------------------
       
        
        %计算KAs的mtau htau
        [tauh_KAs_MSN_Dendrite, taum_KAs_MSN_Dendrite] = tau_KAs_MSN(y(2));
        
        %计算mInf hInf
        mInf_KAs_MSN_Dendrite = infinite(y(2),KAs_MSN_Dendrite.mVhalf,KAs_MSN_Dendrite.mk);
        hInf_KAs_MSN_Dendrite = infinite(y(2),KAs_MSN_Dendrite.hVhalf,KAs_MSN_Dendrite.hk);
        
        %计算 37 38 m,h
        dy(37) = (mInf_KAs_MSN_Dendrite-y(37))/taum_KAs_MSN_Dendrite;
        dy(38) = (hInf_KAs_MSN_Dendrite-y(38))/tauh_KAs_MSN_Dendrite;
        I_KAs_MSN_Dendrite = KAs_MSN_Dendrite.gmax*(y(37)^2)*(KAs_MSN_Dendrite.a*y(38)+(1-KAs_MSN_Dendrite.a))*(y(2)-KAs_MSN_Dendrite.rev);%---------mA
        
        %----------------------KIR----------------------------

        
        %计算KIR的mtau
        taum_KIR_MSN_Dendrite = tau_KIR_MSN(y(2));
        
        %计算mInf
        mInf_KIR_MSN_Dendrite = infinite(y(2),KIR_MSN.mVhalf,KIR_MSN.mk);
        
        %计算 39 m
        dy(39) = (mInf_KIR_MSN_Dendrite-y(39))/taum_KIR_MSN_Dendrite;
        I_KIR_MSN_Dendrite = KIR_MSN.gmax*y(39)*(y(2)-KIR_MSN.rev);%---------mA

        %----------------------CaL1_2----------------------------
        
        %计算CaL1_2的mtau htau
        [tauh_CaL_MSN_Dendrite, taum_CaL_MSN_Dendrite] = tau_CaL_MSN(y(2));
        
        %计算mInf hInf
        mInf_CaL1_2_MSN_Dendrite = infinite(y(2),CaL1_2_MSN.mVhalf,CaL1_2_MSN.mk);
        hInf_CaL1_2_MSN_Dendrite = infinite(y(2),CaL1_2_MSN.mVhalf,CaL1_2_MSN.mk);
        
        %计算 40 41 m,h  % 51 L型Ca CaL_i
        dy(40) = (mInf_CaL1_2_MSN_Dendrite-y(40))/taum_CaL_MSN_Dendrite;
        dy(41) = (hInf_CaL1_2_MSN_Dendrite-y(41))/tauh_CaL_MSN_Dendrite;
        
        g_CaL1_2_MSN_Dendrite = g_Ca_MSN(y(2),y(51),CaL_0_MSN);
        I_CaL1_2_MSN_Dendrite = g_CaL1_2_MSN_Dendrite*CaL1_2_MSN.Pbar*y(40)*y(40)*(y(41)*CaL1_2_MSN.a+(1-CaL1_2_MSN.a));%------mA
        
        %----------------------CaL1_3----------------------------
        
        %计算CaL1_3的mtau htau-----------与CaL1_2的相同
        
        %计算mInf hInf
        mInf_CaL1_3_MSN_Dendrite = infinite(y(2),CaL1_3_MSN.mVhalf,CaL1_3_MSN.mk);
        hInf_CaL1_3_MSN_Dendrite = infinite(y(2),CaL1_3_MSN.mVhalf,CaL1_3_MSN.mk);
        
        %计算 42 43 m,h  % 51 L型Ca CaL_i
        dy(42) = (mInf_CaL1_3_MSN_Dendrite-y(42))/taum_CaL_MSN_Dendrite;
        dy(43) = (hInf_CaL1_3_MSN_Dendrite-y(43))/tauh_CaL_MSN_Dendrite;
        
        g_CaL1_3_MSN_Dendrite = g_Ca_MSN(y(2),y(51),CaL_0_MSN);
        I_CaL1_3_MSN_Dendrite = g_CaL1_3_MSN_Dendrite*CaL1_3_MSN.Pbar*y(42)*y(42)*y(43);%------mA

        
         
        %----------------------CaN----------------------------
        
        
        %计算CaN的mtau htau
        [tauh_CaN_MSN_Dendrite, taum_CaN_MSN_Dendrite] = tau_CaN_MSN(y(2));


        
        %计算mInf hInf
        mInf_CaN_MSN_Dendrite = infinite(y(2), CaN_MSN.mVhalf, CaN_MSN.mk);
        hInf_CaN_MSN_Dendrite = infinite(y(2), CaN_MSN.hVhalf, CaN_MSN.hk);
        
        %计算 44 45 m,h  % 52 Ca Ca_i
        dy(44) = (mInf_CaN_MSN_Dendrite-y(44))/taum_CaN_MSN_Dendrite;
        dy(45) = (hInf_CaN_MSN_Dendrite-y(45))/tauh_CaN_MSN_Dendrite;
        
        g_CaN_MSN_Dendrite = g_Ca_MSN(y(2),y(52),Ca_0_MSN);
        I_CaN_MSN_Dendrite = g_CaN_MSN_Dendrite*CaN_MSN.Pbar*y(44)*y(44)*(CaN_MSN.a*y(45)+(1-CaN_MSN.a));%------mA

             
        %----------------------CaQ----------------------------
        
        
        %计算CaQ的mtau
        taum_CaQ_MSN_Dendrite = tau_CaQ_MSN(y(2));
        
        %计算mInf hInf
        mInf_CaQ_MSN_Dendrite = infinite(y(2),CaQ_MSN.mVhalf,CaQ_MSN.mk);

        
        %计算 46 m  % 52 Ca Ca_i
        dy(46) = (mInf_CaQ_MSN_Dendrite-y(46))/taum_CaQ_MSN_Dendrite;
        
        g_CaQ_MSN_Dendrite = g_Ca_MSN(y(2),y(52),Ca_0_MSN);
        I_CaQ_MSN_Dendrite = g_CaQ_MSN_Dendrite*CaQ_MSN.Pbar*y(46)*y(46);%------mA
        
        
        
        %----------------------CaR----------------------------
        

        
        %计算CaR的mtau htau
        [tauh_CaR_MSN_Dendrite, taum_CaR_MSN_Dendrite] = tau_CaR_MSN(y(2));
        
        %计算mInf hInf
        mInf_CaR_MSN_Dendrite = infinite(y(2),CaR_MSN.mVhalf,CaR_MSN.mk);
        hInf_CaR_MSN_Dendrite = infinite(y(2),CaR_MSN.hVhalf,CaR_MSN.hk);
        
        %计算 47 48 m,h  % 52 Ca Ca_i
        dy(47) = (mInf_CaR_MSN_Dendrite-y(47))/taum_CaR_MSN_Dendrite;
        dy(48) = (hInf_CaR_MSN_Dendrite-y(48))/tauh_CaR_MSN_Dendrite;
        
        g_CaR_MSN_Dendrite = g_Ca_MSN(y(2),y(52),Ca_0_MSN);
        I_CaR_MSN_Dendrite = g_CaR_MSN_Dendrite*CaR_MSN.Pbar*y(47)*y(47)*y(47)*y(48);%------mA
         
        
        %----------------------CaT----------------------------
        
        %计算CaT的mtau htau
        [tauh_CaT_MSN_Dendrite, taum_CaT_MSN_Dendrite] = tau_CaT_MSN(y(2));
        %计算mInf hInf
        mInf_CaT_MSN_Dendrite = infinite(y(2),CaT_MSN.mVhalf,CaT_MSN.mk);
        hInf_CaT_MSN_Dendrite = infinite(y(2),CaT_MSN.hVhalf,CaT_MSN.hk);
        
        %计算 49 50 m,h  % 51 CaL CaL_i
        dy(49) = (mInf_CaT_MSN_Dendrite-y(49))/taum_CaT_MSN_Dendrite;
        dy(50) = (hInf_CaT_MSN_Dendrite-y(50))/tauh_CaT_MSN_Dendrite;
        
        g_CaT_MSN_Dendrite = g_Ca_MSN(y(2),y(51),Ca_0_MSN);
        I_CaT_MSN_Dendrite = g_CaT_MSN_Dendrite*CaT_MSN.Pbar*y(49)*y(49)*y(49)*y(50);%------mA
 
        %----------------------51 CaLdyn ----------------------------        

        I_CaL_MSN_Dendrite = I_CaL1_2_MSN_Dendrite+I_CaL1_3_MSN_Dendrite+I_CaT_MSN_Dendrite;
        drive_channel_CaL_MSN_Dendrite = -drive_MSN * I_CaL_MSN_Dendrite /(2*F_MSN*depth_MSN);%----------把传入的钙（从通道）转换成内部的浓度变化
        %---------------------如果<=0,不能泵入
        if drive_channel_CaL_MSN_Dendrite<=0
            drive_channel_CaL_MSN_Dendrite = 0;
        end
        drive_pump_CaL_MSN_Dendrite = -Kt_Ca_MSN * y(51)/(y(51)+Kd_Ca_MSN);
        dy(51) = drive_channel_CaL_MSN_Dendrite + pump_MSN*drive_pump_CaL_MSN_Dendrite +(CaL_i_inf_MSN-y(51))/tau_R_MSN;
        
        
        %----------------------52 Ca ----------------------------        
        I_Ca_MSN_Dendrite = I_CaN_MSN_Dendrite+I_CaQ_MSN_Dendrite+I_CaR_MSN_Dendrite;
        drive_channel_Ca_MSN_Dendrite = -drive_MSN * I_Ca_MSN_Dendrite /(2*F_MSN*depth_MSN);%----------把传入的钙（从通道）转换成内部的浓度变化
        if drive_channel_Ca_MSN_Dendrite<=0
            drive_channel_Ca_MSN_Dendrite = 0;
        end
        drive_pump_Ca_MSN_Dendrite = -Kt_Ca_MSN * y(52)/(y(52)+Kd_Ca_MSN);
        dy(52) = drive_channel_Ca_MSN_Dendrite + pump_MSN*drive_pump_Ca_MSN_Dendrite +(Ca_i_inf_MSN-y(52))/tau_R_MSN;
        
        %----------------------53 SKKCa ----------------------------
        [tau_SKKCa_MSN_Dendrite,oinf_SKKCa_MSN_Dendrite] = rate_SKKCa_MSN(y(2),y(52));
        dy(53) = (oinf_SKKCa_MSN_Dendrite-y(53))/tau_SKKCa_MSN_Dendrite;
        I_SKKCa_MSN_Dendrite = SKKCa_MSN.gmax*y(53)*(y(2)-SKKCa_MSN.rev);

        
        %----------------------BKKCa ----------------------------
        %---------------54 55 56 :BKKCa    54:open,55:close,56:inactive
        [k1_BKKCa_MSN_Dendrite,k2_BKKCa_MSN_Dendrite,k3_BKKCa_MSN_Dendrite,k4_BKKCa_MSN_Dendrite] = rates_BKKCa_MSN(y(2),y(52));
        dy(54) = -k4_BKKCa_MSN_Dendrite*y(54)-k1_BKKCa_MSN_Dendrite*y(54)+k3_BKKCa_MSN_Dendrite*y(55);%%%%%%%%open state
        dy(55) = -k3_BKKCa_MSN_Dendrite*y(55)+k4_BKKCa_MSN_Dendrite*y(54)+k2_BKKCa_MSN_Dendrite*y(56);%%%%%%%%close state
        dy(56) = k1_BKKCa_MSN_Dendrite*y(54)-k2_BKKCa_MSN_Dendrite*y(56);%%%%%%%inactive state
        I_BKKCa_MSN_Dendrite = BKKCa_MSN.gmax*y(54)*(y(2)-BKKCa_MSN.rev);
        
        
        
        %----------------------I_L漏电流----------------------------             
        I_L_MSN_Dendrite = Leak_MSN.gmax*(y(2)-Leak_MSN.rev);
        
        
        NMDA_Mg_Block_MSN = 1/(1+exp(-(y(1)+15)/16.13));
        T_Presynapse_Pyra2MSN = 1/(1+exp(-(5/5)));
        dy(57) = AMPA_PFCPyra2MSN.alpha*T_Presynapse_Pyra2MSN*...
            (1-y(57))-AMPA_PFCPyra2MSN.beta*y(57);%AMPA Gates
        I_AMPA_Pyra2MSN = AMPA_PFCPyra2MSN.gmax*y(57)*(y(1)-AMPA_PFCPyra2MSN.rev);%Postsynapse:Soma of Pyramidal
        dy(58) = NMDA_PFCPyra2MSN.alpha*T_Presynapse_Pyra2MSN*(1-y(58))-NMDA_PFCPyra2MSN.beta*y(58);%NMDA Gates
        I_NMDA_Pyra2MSN = NMDA_PFCPyra2MSN.gmax*y(58)*NMDA_Mg_Block_MSN*(y(1)-NMDA_PFCPyra2MSN.rev);

        I_MSN_Excitatory = I_AMPA_Pyra2MSN+I_NMDA_Pyra2MSN;
        
        
        %------------主方程-------------
        I_ion_MSN_Soma = (I_NaF_MSN_Soma+I_NaP_MSN_Soma+I_KAf_MSN_Soma+I_KAs_MSN_Soma+I_KIR_MSN_Soma+...
            I_KRP_MSN_Soma+I_CaL1_2_MSN_Soma+I_CaL1_3_MSN_Soma+I_CaN_MSN_Soma+...
            I_CaQ_MSN_Soma+I_CaR_MSN_Soma+I_CaT_MSN_Soma+I_L_MSN_Soma+I_SKKCa_MSN_Soma+I_BKKCa_MSN_Soma)*1000;%1000
        I_ion_MSN_Dendrite = (I_NaF_MSN_Dendrite+I_NaP_MSN_Dendrite+I_KAf_MSN_Dendrite+I_KAs_MSN_Dendrite+I_KIR_MSN_Dendrite+...
            I_CaL1_2_MSN_Dendrite+I_CaL1_3_MSN_Dendrite+I_CaN_MSN_Dendrite+...
            I_CaQ_MSN_Dendrite+I_CaR_MSN_Dendrite+I_CaT_MSN_Dendrite+I_L_MSN_Dendrite+I_SKKCa_MSN_Dendrite+I_BKKCa_MSN_Dendrite)*1000;
        
        dy(1) = -I_ion_MSN_Soma - 10*g_Rall_MSN(1)*(y(1)-y(2)) + I_Ext_MSN_Soma + I_MSN_Excitatory;
        dy(2) = -I_ion_MSN_Dendrite - g_Rall_MSN(2)*(y(2)-y(1)) + I_Ext_MSN_Dendrite;
        

    end



end
