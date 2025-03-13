function [t,y] = NAC_PFC_Dopamine_Ratio_Change(Seed_Net, Seed_Ext, Dopamine_Ratio, I_Ext_MSN_Soma, I_Ext_MSN_Dendrite, I_Ext_PFC_Pyra_Soma,I_Ext_PFC_Pyra_Proximal,I_Ext_PFC_Pyra_Distal, I_Ext_PFC_Inter_CB,I_Ext_PFC_Inter_PV, I_Ext_NAc_Inter_CB, I_Ext_NAc_Inter_PV)
%Dopamine Ratio Change: Dopamine_Ratio ranges from a to b.


% include NAc: 1*MSNs + 1*PV + 1*CB
%         PFC: 20*Pyras + 3*PV + 2*CB
%         In total: 1+1+1+20+3+2=28
% NAc 
% MSN:  1:Soma; 2:Dendrite ;3-56:Ion Channels; 57:PV GABAa Soma;
%       58:CB GABAa Soma; 59-98: 20*(AMPA+NMDA), Soma, PFC
%       In total: 98;
% PV:   99:Soma; 100-101: Ion Channels; 102:CB GABAa Soma; 
%       103-142: 20*Glu(AMPA+NMDA), Soma, PFC;
% CB:   143:Soma, 144-147: Ion Channels; 148:PV GABAa Soma; 
%       149-188: 20*Glu(AMPA+NMDA), Soma, PFC;
%PFC
%Pyra:  189:Soma, 190:Proximal, 191:Distal; 192-200:Ion Channels;
%       201-238: 19*Pyras,AMPA+NMDA, Soma;
%       239-241 3*PV,GABAa Soma; 
%       242-243: 2*CB,GABAa Soma 12+5
%       In total: 20*55=1100, 189-1288: Pyras
%PV:    1289:Soma; 1290-1291: Ion Channels; 1292-1293: 2*PV, GABAa Soma;
%       1294-1295: 2*CB, GABAa Soma; 1296-1335: 20*PV,Glus;
%       In total: 3*47 = 141, 1289-1429,PVs
%CB:    1430:Soma, 1431-1434: Ion Channels, 1435-1437: 3*PV, GABAa Soma;
%       1438: 1*CB, GABAa Soma; 1439-1478: 20*PV,Glus;
%       In total: 2*49 = 98, 1430-1527: CBs


%% 
tic
rng(Seed_Net); %set seed(25);
Net_Strength = rand(28,28); % 1:5,MSN; 6: NAc PV,7: NAc,CB
                            % 8-27: Pyras; 28-30: PFC PVs; 31-32: PFC CBs;
                            % a(i,j) = a(post,pre)
Net_Strength = Net_Strength-diag(diag(Net_Strength));%diag=0
Net_Strength = roundn(Net_Strength,-4);

% Parameters
% MSN
V_MSN_Soma = -90;
V_MSN_Dendrite = -90;
V_PFC_Inter_PV = -64;
V_PFC_Inter_CB = -64;
V_PFC_Pyra_Soma = -64.8;
V_PFC_Pyra_Proximal = -64;
V_PFC_Pyra_Distal = -64;
y0 = Initial_NAc_PFC(V_PFC_Inter_PV,V_PFC_Pyra_Soma,V_PFC_Pyra_Proximal,V_PFC_Pyra_Distal,...
    V_MSN_Soma,V_MSN_Dendrite,V_PFC_Inter_CB, Dopamine_Ratio(1));



%% Parameters
NaF_MSN_Soma = MSN_ionmodel1('NaF', 1.5, nan, -23.9, -62.9, -11.8, 10.7, 50);
NaF_MSN_Dendrite = MSN_ionmodel1('NaF', 0.0195, nan, -23.9, -62.9, -11.8, 10.7, 50);
NaP_MSN_Soma = MSN_ionmodel1('NaP', 4e-5, nan, -52.6, -48.8, -4.6, 10, 50);
NaP_MSN_Dendrite = MSN_ionmodel1('NaP', 1.38e-7, nan, -52.6, -48.8, -4.6, 10, 50);
KAf_MSN_Soma = MSN_ionmodel1('KAf', 0.225, nan, -10, -75.6, -17.7, 10, -90);
% KAf_MSN_Soma = MSN_ionmodel1('KAf', 0.3, nan, -10, -75.6, -17.7, 10, -90);

KAf_MSN_Dendrite = MSN_ionmodel1('KAf', 0.021, nan, -10, -75.6, -17.7, 10, -90);
% KAs_MSN_Soma = MSN_ionmodel1('KAs', 0.0104, 0.996, -27, -33.5, -16, 21.5, -90);% Standard
% KAs_MSN_Soma = MSN_ionmodel1('KAs', 0.0104/2, 0.996, -27, -33.5, -16, 21.5, -90);% Dopamine
KAs_MSN_Soma = MSN_ionmodel1('KAs', (0.0104-0.0104/2*Dopamine_Ratio), 0.996, -27, -33.5, -16, 21.5, -90);% Dopamine Ratio
% KAs_MSN_Dendrite = MSN_ionmodel1('KAs', 9.51e-4, 0.996, -27, -33.5, -16, 21.5, -90);% Standard
% KAs_MSN_Dendrite = MSN_ionmodel1('KAs', 9.51e-4/2, 0.996, -27, -33.5, -16, 21.5, -90);% Dopamine
KAs_MSN_Dendrite = MSN_ionmodel1('KAs', (9.51e-4-9.51e-4/2*Dopamine_Ratio), 0.996, -27, -33.5, -16, 21.5, -90);% Dopamine Ratio

KIR_MSN = MSN_ionmodel1('KIR', 1.4e-4, nan, -82, nan, 13, nan, -90);
% KRP_MSN = MSN_ionmodel1('KRP', 0.001, 0.7, -13.5, -54.7, -11.8, 18.6, -90);
KRP_MSN = MSN_ionmodel1('KRP', 0.001, 0.7, -13.5, -54.7, -11.8, 18.6, -90);
Leak_MSN = MSN_ionmodel1('Leak', 11.5e-6, nan, nan, nan, nan, nan, -70);
BKKCa_MSN = MSN_ionmodel1('BKKCa', 0.001, nan, nan, nan, nan, nan, -90);
SKKCa_MSN = MSN_ionmodel1('SKKCa', 0.145, nan, nan, nan, nan, nan, -90);
% CaL1_2_MSN = MSN_ionmodel2('CaL1.2', 6.7e-6, 0.17, -8.9, -13.4, -6.7, 11.9);%Standard
% CaL1_2_MSN = MSN_ionmodel2('CaL1.2', 6.7e-6*0.8, 0.17, -8.9, -13.4, -6.7, 11.9);%Dopamine
CaL1_2_MSN = MSN_ionmodel2('CaL1.2', (6.7e-6-6.7e-6*0.2*Dopamine_Ratio), 0.17, -8.9, -13.4, -6.7, 11.9);%Dopamine Ratio

CaL1_3_MSN = MSN_ionmodel2('CaL1.3', 4.25e-7, nan, -33, -13.4, -6.7, 11.9);
CaN_MSN = MSN_ionmodel2('CaN', 1.0e-5, 0.21, -8.7, -74.8, -7.4, 6.5);
CaQ_MSN = MSN_ionmodel2('CaQ', 6.0e-6, nan, -9.0, nan, -6.6, nan);
CaR_MSN = MSN_ionmodel2('CaR', 2.6e-5, nan, -10.3, -33.3, -6.6, 17);
CaT_MSN = MSN_ionmodel2('CaT', 4e-7, nan, -51.73, -80, -6.53, 6.7);
g_Rall_MSN = MSN_RallModel;

Na_PFC_Inter_PV = PFC_Inter_ionmodel('Na', 35, 5, 55);
K_PFC_Inter_PV = PFC_Inter_ionmodel('K', 9, 5, -90);
Leak_PFC_Inter_PV = PFC_Inter_ionmodel('Leak', 0.1, 1, -65);
Na_PFC_Inter_CB = PFC_Inter_ionmodel('Na', 35, 5, 55);
K_PFC_Inter_CB = PFC_Inter_ionmodel('K', 9, 5, -85);
Ca_PFC_Inter_CB = PFC_Pyra_ionmodel('Ca', 1, 1, 120);
KCa_PFC_Inter_CB = PFC_Pyra_ionmodel('KCa', 1, 1, -85);
h_PFC_Inter_CB = PFC_Pyra_ionmodel('KCa', 0.15, 1, -40);
Leak_PFC_Inter_CB = PFC_Pyra_ionmodel('Leak', 0.1, 1, -65);

Na_PFC_Pyra_Soma = PFC_Pyra_ionmodel('Na', 55, 4, 55);
K_PFC_Pyra_Soma = PFC_Pyra_ionmodel('K', 60, 4, -80);
% Ca_PFC_Pyra_Soma = PFC_Pyra_ionmodel('Ca', 1.5, 1, 120); % Standard
% Ca_PFC_Pyra_Soma = PFC_Pyra_ionmodel('Ca', 1.2, 1, 120);% Dopamine
Ca_PFC_Pyra_Soma = PFC_Pyra_ionmodel('Ca', (1.5-0.3*Dopamine_Ratio), 1, 120);% Dopamine Ratio

% CaN_PFC_Pyra_Soma = PFC_Pyra_ionmodel('CaN', 0.025, 1, -20);% Standard
% CaN_PFC_Pyra_Soma = PFC_Pyra_ionmodel('CaN', 0.02, 1, -20);% Dopamine
CaN_PFC_Pyra_Soma = PFC_Pyra_ionmodel('CaN', (0.025-0.005*Dopamine_Ratio), 1, -20);% Dopamine Ratio


Leak_PFC_Pyra = PFC_Pyra_ionmodel('Leak', 0.05, 1, -70);
NaP_PFC_Pyra_Proximal = PFC_Pyra_ionmodel('NaP', 0.15, 1, 55);
% KS_PFC_Pyra_Proximal = PFC_Pyra_ionmodel('KS', 16, 1, -80); % g = 16 burst, standard=2
% KS_PFC_Pyra_Proximal = PFC_Pyra_ionmodel('KS', 8, 1, -80); %Dopamine
KS_PFC_Pyra_Proximal = PFC_Pyra_ionmodel('KS', (16-8*Dopamine_Ratio), 1, -80); %Dopamine Ratio


A_PFC_Pyra_Distal = PFC_Pyra_ionmodel('A', 1, 1, -80);
% Ca_PFC_Pyra_Distal = PFC_Pyra_ionmodel('Ca', 0.25, 1, 120);% Standard
% Ca_PFC_Pyra_Distal = PFC_Pyra_ionmodel('Ca', 0.2, 1, 120); % Dopamine
Ca_PFC_Pyra_Distal = PFC_Pyra_ionmodel('Ca', (0.25-0.05*Dopamine_Ratio), 1, 120); % Dopamine Ratio





%Synapse
AMPA_PFCPyra2PFCPV = SynapseChannel_Debug('AMPA','PFCPyra','PFCPV',0.0003125*10,1,0.2,0);
NMDA_PFCPyra2PFCPV = SynapseChannel_Debug('NMDA','PFCPyra','PFCPV',0.0003125*10,0.072,0.0067,0);
GABAa_PFCInter2PFCPV = SynapseChannel_Debug('GABAa','PFCInter','PFCPV',0.0025*10*3,5,0.18,-80);

AMPA_PFCPyra2PFCCB = SynapseChannel_Debug('AMPA','PFCPyra','PFCCB',0.0003125*1,1,0.2,0);
NMDA_PFCPyra2PFCCB = SynapseChannel_Debug('NMDA','PFCPyra','PFCCB',0.0003125*1,0.072,0.0067,0);
GABAa_PFCInter2PFCCB = SynapseChannel_Debug('GABAa','PFCInter','PFCCB',0.0025*10*3,5,0.18,-80);

AMPA_PFCPyra2PFCPyra = SynapseChannel_Debug('AMPA','PFCPyra','PFCPyra',0.0013*10,1,0.2,0);
NMDA_PFCPyra2PFCPyra = SynapseChannel_Debug('NMDA','PFCPyra','PFCPyra',0.0013*10,0.072,0.0067,0);
GABAa_PFCInter2PFCPyra = SynapseChannel_Debug('GABAa','PFCInter','PFCPyra',0.005*10*2,5,0.18,-80);

AMPA_PFCPyra2MSN = SynapseChannel_Debug('AMPA','PFCPyra','MSN',0.002*10,1,0.2,0);
NMDA_PFCPyra2MSN = SynapseChannel_Debug('NMDA','PFCPyra','MSN',0.002*10,0.072,0.0067,0);
GABAa_NAcInter2MSN = SynapseChannel_Debug('GABAa','NAcInter','MSN',0.01*10*3,5,0.18,-80);
% AMPA_PFCPyra2MSN = SynapseChannel('AMPA','PFCPyra','MSN',0.002,1,0.2,0);
% NMDA_PFCPyra2MSN = SynapseChannel('NMDA','PFCPyra','MSN',0.002,0.072,0.0067,0);
% GABAa_NAcInter2MSN = SynapseChannel('GABAa','NAcInter','MSN',0.01,5,0.18,-80);

AMPA_PFCPyra2NAcPV = SynapseChannel_Debug('AMPA','PFCPyra','NAcPV',0.0003125*10,1,0.2,0);
NMDA_PFCPyra2NAcPV = SynapseChannel_Debug('NMDA','PFCPyra','NAcPV',0.0003125*10,0.072,0.0067,0);
GABAa_NAcInter2NAcPV = SynapseChannel_Debug('GABAa','NAcInter','NAcPV',0.01*10*3,5,0.18,-80);

AMPA_PFCPyra2NAcCB = SynapseChannel_Debug('AMPA','PFCPyra','NAcCB',0.0003125*1,1,0.2,0);
NMDA_PFCPyra2NAcCB = SynapseChannel_Debug('NMDA','PFCPyra','NAcCB',0.0003125*1,0.072,0.0067,0);
GABAa_NAcInter2NAcCB = SynapseChannel_Debug('GABAa','NAcInter','NAcCB',0.001*10*3,5,0.18,-80);


%% 
t = 0:0.02:3000;
t0 = 0:0.02:3005;
rng(Seed_Ext)
% [t,y] = NAC_PFC_Change_Ext(5,4,2,2,2,0.8,0.8,0.8,1);
I_Ext_MSN_Soma_Real = roundn(I_Ext_MSN_Soma*ones(size(t0)) + 2*rand(size(t0)),-4);
I_Ext_MSN_Dendrite_Real = roundn(I_Ext_MSN_Dendrite*ones(size(t0)) + 2*rand(size(t0)),-4);
I_Ext_PFC_Pyra_Soma_Real = roundn(I_Ext_PFC_Pyra_Soma*ones(size(t0)) + 1*rand(size(t0)),-4);
I_Ext_PFC_Pyra_Proximal_Real = roundn(I_Ext_PFC_Pyra_Proximal*ones(size(t0)) + 1*rand(size(t0)),-4);
I_Ext_PFC_Pyra_Distal_Real = roundn(I_Ext_PFC_Pyra_Distal*ones(size(t0)) + 1*rand(size(t0)),-4);
I_Ext_PFC_Inter_CB_Real = roundn(I_Ext_PFC_Inter_CB*ones(size(t0)) + 0.4*rand(size(t0)),-4);
I_Ext_PFC_Inter_PV_Real = roundn(I_Ext_PFC_Inter_PV*ones(size(t0)) + 0.4*rand(size(t0)),-4);
I_Ext_NAc_Inter_CB_Real = roundn(I_Ext_NAc_Inter_CB*ones(size(t0)) + 0.4*rand(size(t0)),-4);
I_Ext_NAc_Inter_PV_Real = roundn(I_Ext_NAc_Inter_PV*ones(size(t0)) + 1*rand(size(t0)),-4);

[t,y] = ode23(@S_NAC_PFC_Change_Ext,t,y0);
toc
%%
    function dy = S_NAC_PFC_Change_Ext(t,y)

        %%
        dy = zeros(1527,1);
        %% Potential
        % NAc MSN:1 , PV: 99 CB:143
        % PFC Pyras:189+55*0~19; PV:1289+47*0~2; CB:1430+49*0~1

       %% NAc, MSN
       % NAc
       % MSN:  1:Soma; 2:Dendrite ;3-56:Ion Channels; 57:PV GABAa Soma;
       %       58:CB GABAa Soma; 59-98: 20*(AMPA+NMDA), Soma, PFC
       %       In total: 98;
           
       %----------------------NaF----------------------------
       
       % NaF, mtau htau
       [tauh_NaF_MSN_Soma, taum_NaF_MSN_Soma] = tau_NaF_MSN(y(1));
       
       % mInf hInf
       mInf_NaF_MSN_Soma = infinite(y(1),NaF_MSN_Soma.mVhalf,NaF_MSN_Soma.mk);
       hInf_NaF_MSN_Soma = infinite(y(1),NaF_MSN_Soma.hVhalf,NaF_MSN_Soma.hk);
       
       % 3 4 m,h
       dy(3) = (mInf_NaF_MSN_Soma-y(3))/taum_NaF_MSN_Soma;
       dy(4) = (hInf_NaF_MSN_Soma-y(4))/tauh_NaF_MSN_Soma;
       I_NaF_MSN_Soma = NaF_MSN_Soma.gmax*(y(3)^3)*y(4)*(y(1)-NaF_MSN_Soma.rev);%-----mA
       
       %----------------------NaP----------------------------
       
       % NaP, mtau htau
       [tauh_NaP_MSN_Soma, taum_NaP_MSN_Soma] = tau_NaP_MSN(y(1));
       
       % mInf hInf
       mInf_NaP_MSN_Soma = infinite(y(1),NaP_MSN_Soma.mVhalf,NaP_MSN_Soma.mk);
       hInf_NaP_MSN_Soma = infinite(y(1),NaP_MSN_Soma.hVhalf,NaP_MSN_Soma.hk);
       
       % 5 6  m,h
       dy(5) = (mInf_NaP_MSN_Soma-y(5))/taum_NaP_MSN_Soma;
       dy(6) = (hInf_NaP_MSN_Soma-y(6))/tauh_NaP_MSN_Soma;
       I_NaP_MSN_Soma = NaP_MSN_Soma.gmax*y(5)*y(6)*(y(1)-NaP_MSN_Soma.rev);%---------mA
       
       %----------------------KAf----------------------------
       
       % KAf, mtau htau
       [tauh_KAf_MSN_Soma, taum_KAf_MSN_Soma] = tau_KAf_MSN(y(1));
       
       % mInf hInf
       mInf_KAf_MSN_Soma = infinite(y(1),KAf_MSN_Soma.mVhalf,KAf_MSN_Soma.mk);
       hInf_KAf_MSN_Soma = infinite(y(1),KAf_MSN_Soma.hVhalf,KAf_MSN_Soma.hk);
       
       % 7 8 m,h
       dy(7) = (mInf_KAf_MSN_Soma-y(7))/taum_KAf_MSN_Soma;
       dy(8) = (hInf_KAf_MSN_Soma-y(8))/tauh_KAf_MSN_Soma;
       I_KAf_MSN_Soma = KAf_MSN_Soma.gmax*(y(7)^2)*y(8)*(y(1)-KAf_MSN_Soma.rev);%---------mA
       
       
       %----------------------KAs----------------------------
       
       
       % KAs, mtau htau
       [tauh_KAs_MSN_Soma, taum_KAs_MSN_Soma] = tau_KAs_MSN(y(1));
       
       % mInf hInf
       mInf_KAs_MSN_Soma = infinite(y(1),KAs_MSN_Soma.mVhalf,KAs_MSN_Soma.mk);
       hInf_KAs_MSN_Soma = infinite(y(1),KAs_MSN_Soma.hVhalf,KAs_MSN_Soma.hk);
       
       % 9 10 m,h
       dy(9) = (mInf_KAs_MSN_Soma-y(9))/taum_KAs_MSN_Soma;
       dy(10) = (hInf_KAs_MSN_Soma-y(10))/tauh_KAs_MSN_Soma;
       I_KAs_MSN_Soma = KAs_MSN_Soma.gmax(floor(t/0.02)+1)*(y(9)^2)*(KAs_MSN_Soma.a*y(10)+(1-KAs_MSN_Soma.a))*(y(1)-KAs_MSN_Soma.rev);%---------mA
       
       %----------------------KIR----------------------------
       
       
       % KIR, mtau
       taum_KIR_MSN_Soma = tau_KIR_MSN(y(1));
       
       % mInf
       mInf_KIR_MSN_Soma = infinite(y(1),KIR_MSN.mVhalf,KIR_MSN.mk);
       
       % 11 m
       dy(11) = (mInf_KIR_MSN_Soma-y(11))/taum_KIR_MSN_Soma;
       I_KIR_MSN_Soma = KIR_MSN.gmax*y(11)*(y(1)-KIR_MSN.rev);%---------mA
       
       %----------------------KRP----------------------------
       
       % KRP, mtau htau
       [tauh_KRP_MSN_Soma, taum_KRP_MSN_Soma] = tau_KRP_MSN(y(1));
       
       % mInf hInf
       mInf_KRP_MSN_Soma = infinite(y(1),KRP_MSN.mVhalf,KRP_MSN.mk);
       hInf_KRP_MSN_Soma = infinite(y(1),KRP_MSN.mVhalf,KRP_MSN.mk);
       
       % 12 13 m,h
       dy(12) = (mInf_KRP_MSN_Soma-y(12))/taum_KRP_MSN_Soma;
       dy(13) = (hInf_KRP_MSN_Soma-y(13))/tauh_KRP_MSN_Soma;
       I_KRP_MSN_Soma = KRP_MSN.gmax*y(12)*(KRP_MSN.a*y(13)+(1-KRP_MSN.a))*(y(1)-KRP_MSN.rev);%---------mA
       
       %----------------------Ca dynamics parameters----------------------------
       
       z_MSN = 2;
       F_MSN = 96489;%-------------Faraday(C/mod)
       R_MSN = 8.31;%------------J/(mol*K)
       T_MSN = 35+273.15;%--------K
       Ca_0_MSN = 5;%---------Ca store
       CaL_0_MSN = 5;%--------CaL store
       
       
       
       
       %----------------------CaL1_2----------------------------
       
       % CaL1_2, mtau htau
       [tauh_CaL_MSN_Soma, taum_CaL_MSN_Soma] = tau_CaL_MSN(y(1));
       
       % mInf hInf
       mInf_CaL1_2_MSN_Soma = infinite(y(1),CaL1_2_MSN.mVhalf,CaL1_2_MSN.mk);
       hInf_CaL1_2_MSN_Soma = infinite(y(1),CaL1_2_MSN.mVhalf,CaL1_2_MSN.mk);
       
       %  14 15 m,h
       dy(14) = (mInf_CaL1_2_MSN_Soma-y(14))/taum_CaL_MSN_Soma;
       dy(15) = (hInf_CaL1_2_MSN_Soma-y(15))/tauh_CaL_MSN_Soma;
       
       g_CaL1_2_MSN_Soma = g_Ca_MSN(y(1),y(25),CaL_0_MSN);
       I_CaL1_2_MSN_Soma = g_CaL1_2_MSN_Soma*CaL1_2_MSN.Pbar(floor(t/0.02)+1)*y(14)*y(14)*(y(15)*CaL1_2_MSN.a+(1-CaL1_2_MSN.a));%------mA
       
       %----------------------CaL1_3----------------------------
       
       % CaL1_3, mtau htau
       
       % mInf hInf
       mInf_CaL1_3_MSN_Soma = infinite(y(1),CaL1_3_MSN.mVhalf,CaL1_3_MSN.mk);
       hInf_CaL1_3_MSN_Soma = infinite(y(1),CaL1_3_MSN.mVhalf,CaL1_3_MSN.mk);
       
       % 16 17m,h
       dy(16) = (mInf_CaL1_3_MSN_Soma-y(16))/taum_CaL_MSN_Soma;
       dy(17) = (hInf_CaL1_3_MSN_Soma-y(17))/tauh_CaL_MSN_Soma;
       
       g_CaL1_3_MSN_Soma = g_Ca_MSN(y(1),y(25),CaL_0_MSN);
       I_CaL1_3_MSN_Soma = g_CaL1_3_MSN_Soma*CaL1_3_MSN.Pbar*y(16)*y(16)*y(17);%------mA
       
       
       
       %----------------------CaN----------------------------
       
       
       % CaN, mtau htau
       [tauh_CaN_MSN_Soma, taum_CaN_MSN_Soma] = tau_CaN_MSN(y(1));
       
       
       
       % mInf hInf
       mInf_CaN_MSN_Soma = infinite(y(1), CaN_MSN.mVhalf, CaN_MSN.mk);
       hInf_CaN_MSN_Soma = infinite(y(1), CaN_MSN.hVhalf, CaN_MSN.hk);
       
       %  18 19 m,h
       dy(18) = (mInf_CaN_MSN_Soma-y(18))/taum_CaN_MSN_Soma;
       dy(19) = (hInf_CaN_MSN_Soma-y(19))/tauh_CaN_MSN_Soma;
       
       g_CaN_MSN_Soma = g_Ca_MSN(y(1),y(26),Ca_0_MSN);
       I_CaN_MSN_Soma = g_CaN_MSN_Soma*CaN_MSN.Pbar*y(18)*y(18)*(CaN_MSN.a*y(19)+(1-CaN_MSN.a));%------mA
       
       
       %----------------------CaQ----------------------------
       
       
       % CaQ, mtau
       taum_CaQ_MSN_Soma = tau_CaQ_MSN(y(1));
       
       % mInf hInf
       mInf_CaQ_MSN_Soma = infinite(y(1),CaQ_MSN.mVhalf,CaQ_MSN.mk);
       
       
       %  20 m
       dy(20) = (mInf_CaQ_MSN_Soma-y(20))/taum_CaQ_MSN_Soma;
       
       g_CaQ_MSN_Soma = g_Ca_MSN(y(1),y(26),Ca_0_MSN);
       I_CaQ_MSN_Soma = g_CaQ_MSN_Soma*CaQ_MSN.Pbar*y(20)*y(20);%------mA
       
       
       
       %----------------------CaR----------------------------
       
       
       
       % CaR, mtau htau
       [tauh_CaR_MSN_Soma, taum_CaR_MSN_Soma] = tau_CaR_MSN(y(1));
       
       % mInf hInf
       mInf_CaR_MSN_Soma = infinite(y(1),CaR_MSN.mVhalf,CaR_MSN.mk);
       hInf_CaR_MSN_Soma = infinite(y(1),CaR_MSN.hVhalf,CaR_MSN.hk);
       
       %  21 22 m,h
       dy(21) = (mInf_CaR_MSN_Soma-y(21))/taum_CaR_MSN_Soma;
       dy(22) = (hInf_CaR_MSN_Soma-y(22))/tauh_CaR_MSN_Soma;
       
       g_CaR_MSN_Soma = g_Ca_MSN(y(1),y(26),Ca_0_MSN);
       I_CaR_MSN_Soma = g_CaR_MSN_Soma*CaR_MSN.Pbar*y(21)*y(21)*y(21)*y(22);%------mA
       
       
       %----------------------CaT----------------------------
       
       % CaT, mtau htau
       [tauh_CaT_MSN_Soma, taum_CaT_MSN_Soma] = tau_CaT_MSN(y(1));
       % mInf hInf
       mInf_CaT_MSN_Soma = infinite(y(1),CaT_MSN.mVhalf,CaT_MSN.mk);
       hInf_CaT_MSN_Soma = infinite(y(1),CaT_MSN.hVhalf,CaT_MSN.hk);
       
       %  23 24 m,h
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
       drive_channel_CaL_MSN_Soma = -drive_MSN * I_CaL_MSN_Soma /(2*F_MSN*depth_MSN);

       if drive_channel_CaL_MSN_Soma<=0
           drive_channel_CaL_MSN_Soma = 0;
       end
       drive_pump_CaL_MSN_Soma = -Kt_Ca_MSN * y(25)/(y(25)+Kd_Ca_MSN);
       dy(25) = drive_channel_CaL_MSN_Soma + pump_MSN*drive_pump_CaL_MSN_Soma +(CaL_i_inf_MSN-y(25))/tau_R_MSN;
       
       
       %----------------------26 Ca ----------------------------
       I_Ca_MSN_Soma = I_CaN_MSN_Soma+I_CaQ_MSN_Soma+I_CaR_MSN_Soma;
       drive_channel_Ca_MSN_Soma = -drive_MSN * I_Ca_MSN_Soma /(2*F_MSN*depth_MSN);
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
       
       
       
       %----------------------I_L----------------------------
       I_L_MSN_Soma = Leak_MSN.gmax*(y(1)-Leak_MSN.rev);
       
       %----------------------Proximal 2,31-56----------------------------
       %---------------           31,32 NaF m h ; 33,34 NaP m h ; 35,36 KAf m h;
       %---------------           37,38 KAs m h ; 39 KIR m ; 40,41 CaL1_2 m h ;
       %---------------           42,43 CaL1_3 m h ; 44,45 CaN m h ; 46 CaQ m ; 47,48 CaR m h ;
       %---------------           59,50 CaT m h;
       %---------------           51 LÐÍCa CaL_i  ; 52:Ca_i
       %---------------           53 SKKCa
       %---------------           54,55,56:BKKCa    54:open,55:close,56:inactive
       
       %----------------------NaF---------------------------
       
       % NaF, mtau htau
       [tauh_NaF_MSN_Dendrite, taum_NaF_MSN_Dendrite] = tau_NaF_MSN(y(2));
       
       % mInf hInf
       mInf_NaF_MSN_Dendrite = infinite(y(2),NaF_MSN_Dendrite.mVhalf,NaF_MSN_Dendrite.mk);
       hInf_NaF_MSN_Dendrite = infinite(y(2),NaF_MSN_Dendrite.hVhalf,NaF_MSN_Dendrite.hk);
       
       % 31 32 m,h
       dy(31) = (mInf_NaF_MSN_Dendrite-y(31))/taum_NaF_MSN_Dendrite;
       dy(32) = (hInf_NaF_MSN_Dendrite-y(32))/tauh_NaF_MSN_Dendrite;
       I_NaF_MSN_Dendrite = NaF_MSN_Dendrite.gmax*(y(31)^3)*y(32)*(y(2)-NaF_MSN_Dendrite.rev);%-----mA
       
       %----------------------NaP----------------------------
       
       % NaP, mtau htau
       [tauh_NaP_MSN_Dendrite, taum_NaP_MSN_Dendrite] = tau_NaP_MSN(y(2));
       
       % mInf hInf
       mInf_NaP_MSN_Dendrite = infinite(y(2),NaP_MSN_Dendrite.mVhalf,NaP_MSN_Dendrite.mk);
       hInf_NaP_MSN_Dendrite = infinite(y(2),NaP_MSN_Dendrite.hVhalf,NaP_MSN_Dendrite.hk);
       
       %  33 34  m,h
       dy(33) = (mInf_NaP_MSN_Dendrite-y(33))/taum_NaP_MSN_Dendrite;
       dy(34) = (hInf_NaP_MSN_Dendrite-y(34))/tauh_NaP_MSN_Dendrite;
       I_NaP_MSN_Dendrite = NaP_MSN_Dendrite.gmax*y(33)*y(34)*(y(2)-NaP_MSN_Dendrite.rev);%---------mA
       
       %----------------------KAf----------------------------
       
       % KAf, mtau htau
       [tauh_KAf_MSN_Dendrite, taum_KAf_MSN_Dendrite] = tau_KAf_MSN(y(2));
       
       % mInf hInf
       mInf_KAf_MSN_Dendrite = infinite(y(2),KAf_MSN_Dendrite.mVhalf,KAf_MSN_Dendrite.mk);
       hInf_KAf_MSN_Dendrite = infinite(y(2),KAf_MSN_Dendrite.hVhalf,KAf_MSN_Dendrite.hk);
       
       % 35 36 m,h
       dy(35) = (mInf_KAf_MSN_Dendrite-y(35))/taum_KAf_MSN_Dendrite;
       dy(36) = (hInf_KAf_MSN_Dendrite-y(36))/tauh_KAf_MSN_Dendrite;
       I_KAf_MSN_Dendrite = KAf_MSN_Dendrite.gmax*(y(35)^2)*y(36)*(y(2)-KAf_MSN_Dendrite.rev);%---------mA
       
       
       %----------------------KAs----------------------------
       
       
       % KAs, mtau htau
       [tauh_KAs_MSN_Dendrite, taum_KAs_MSN_Dendrite] = tau_KAs_MSN(y(2));
       
       % mInf hInf
       mInf_KAs_MSN_Dendrite = infinite(y(2),KAs_MSN_Dendrite.mVhalf,KAs_MSN_Dendrite.mk);
       hInf_KAs_MSN_Dendrite = infinite(y(2),KAs_MSN_Dendrite.hVhalf,KAs_MSN_Dendrite.hk);
       
       %  37 38 m,h
       dy(37) = (mInf_KAs_MSN_Dendrite-y(37))/taum_KAs_MSN_Dendrite;
       dy(38) = (hInf_KAs_MSN_Dendrite-y(38))/tauh_KAs_MSN_Dendrite;
       I_KAs_MSN_Dendrite = KAs_MSN_Dendrite.gmax(floor(t/0.02)+1)*(y(37)^2)*(KAs_MSN_Dendrite.a*y(38)+(1-KAs_MSN_Dendrite.a))*(y(2)-KAs_MSN_Dendrite.rev);%---------mA
       
       %----------------------KIR----------------------------
       
       
       % KIR, mtau
       taum_KIR_MSN_Dendrite = tau_KIR_MSN(y(2));
       
       % mInf
       mInf_KIR_MSN_Dendrite = infinite(y(2),KIR_MSN.mVhalf,KIR_MSN.mk);
       
       %  39 m
       dy(39) = (mInf_KIR_MSN_Dendrite-y(39))/taum_KIR_MSN_Dendrite;
       I_KIR_MSN_Dendrite = KIR_MSN.gmax*y(39)*(y(2)-KIR_MSN.rev);%---------mA
       
       %----------------------CaL1_2----------------------------
       
       % CaL1_2, mtau htau
       [tauh_CaL_MSN_Dendrite, taum_CaL_MSN_Dendrite] = tau_CaL_MSN(y(2));
       
       % mInf hInf
       mInf_CaL1_2_MSN_Dendrite = infinite(y(2),CaL1_2_MSN.mVhalf,CaL1_2_MSN.mk);
       hInf_CaL1_2_MSN_Dendrite = infinite(y(2),CaL1_2_MSN.mVhalf,CaL1_2_MSN.mk);
       
       %  40 41 m,h
       dy(40) = (mInf_CaL1_2_MSN_Dendrite-y(40))/taum_CaL_MSN_Dendrite;
       dy(41) = (hInf_CaL1_2_MSN_Dendrite-y(41))/tauh_CaL_MSN_Dendrite;
       
       g_CaL1_2_MSN_Dendrite = g_Ca_MSN(y(2),y(51),CaL_0_MSN);
       I_CaL1_2_MSN_Dendrite = g_CaL1_2_MSN_Dendrite*CaL1_2_MSN.Pbar(floor(t/0.02)+1)*y(40)*y(40)*(y(41)*CaL1_2_MSN.a+(1-CaL1_2_MSN.a));%------mA
       
       %----------------------CaL1_3----------------------------
       
       % CaL1_3, mtau htau
       
       % mInf hInf
       mInf_CaL1_3_MSN_Dendrite = infinite(y(2),CaL1_3_MSN.mVhalf,CaL1_3_MSN.mk);
       hInf_CaL1_3_MSN_Dendrite = infinite(y(2),CaL1_3_MSN.mVhalf,CaL1_3_MSN.mk);
       
       %  42 43 m,h
       dy(42) = (mInf_CaL1_3_MSN_Dendrite-y(42))/taum_CaL_MSN_Dendrite;
       dy(43) = (hInf_CaL1_3_MSN_Dendrite-y(43))/tauh_CaL_MSN_Dendrite;
       
       g_CaL1_3_MSN_Dendrite = g_Ca_MSN(y(2),y(51),CaL_0_MSN);
       I_CaL1_3_MSN_Dendrite = g_CaL1_3_MSN_Dendrite*CaL1_3_MSN.Pbar*y(42)*y(42)*y(43);%------mA
       
       
       
       %----------------------CaN----------------------------
       
       
       % CaN, mtau htau
       [tauh_CaN_MSN_Dendrite, taum_CaN_MSN_Dendrite] = tau_CaN_MSN(y(2));
       
       
       
       % mInf hInf
       mInf_CaN_MSN_Dendrite = infinite(y(2), CaN_MSN.mVhalf, CaN_MSN.mk);
       hInf_CaN_MSN_Dendrite = infinite(y(2), CaN_MSN.hVhalf, CaN_MSN.hk);
       
       %  44 45 m,h
       dy(44) = (mInf_CaN_MSN_Dendrite-y(44))/taum_CaN_MSN_Dendrite;
       dy(45) = (hInf_CaN_MSN_Dendrite-y(45))/tauh_CaN_MSN_Dendrite;
       
       g_CaN_MSN_Dendrite = g_Ca_MSN(y(2),y(52),Ca_0_MSN);
       I_CaN_MSN_Dendrite = g_CaN_MSN_Dendrite*CaN_MSN.Pbar*y(44)*y(44)*(CaN_MSN.a*y(45)+(1-CaN_MSN.a));%------mA
       
       
       %----------------------CaQ----------------------------
       
       
       % CaQ, mtau
       taum_CaQ_MSN_Dendrite = tau_CaQ_MSN(y(2));
       
       % mInf hInf
       mInf_CaQ_MSN_Dendrite = infinite(y(2),CaQ_MSN.mVhalf,CaQ_MSN.mk);
       
       
       %  46 m
       dy(46) = (mInf_CaQ_MSN_Dendrite-y(46))/taum_CaQ_MSN_Dendrite;
       
       g_CaQ_MSN_Dendrite = g_Ca_MSN(y(2),y(52),Ca_0_MSN);
       I_CaQ_MSN_Dendrite = g_CaQ_MSN_Dendrite*CaQ_MSN.Pbar*y(46)*y(46);%------mA
       
       
       
       %----------------------CaR----------------------------
       
       
       
       % CaR, mtau htau
       [tauh_CaR_MSN_Dendrite, taum_CaR_MSN_Dendrite] = tau_CaR_MSN(y(2));
       
       % mInf hInf
       mInf_CaR_MSN_Dendrite = infinite(y(2),CaR_MSN.mVhalf,CaR_MSN.mk);
       hInf_CaR_MSN_Dendrite = infinite(y(2),CaR_MSN.hVhalf,CaR_MSN.hk);
       
       %  47 48 m,h
       dy(47) = (mInf_CaR_MSN_Dendrite-y(47))/taum_CaR_MSN_Dendrite;
       dy(48) = (hInf_CaR_MSN_Dendrite-y(48))/tauh_CaR_MSN_Dendrite;
       
       g_CaR_MSN_Dendrite = g_Ca_MSN(y(2),y(52),Ca_0_MSN);
       I_CaR_MSN_Dendrite = g_CaR_MSN_Dendrite*CaR_MSN.Pbar*y(47)*y(47)*y(47)*y(48);%------mA
       
       
       %----------------------CaT----------------------------
       
       % CaT, mtau htau
       [tauh_CaT_MSN_Dendrite, taum_CaT_MSN_Dendrite] = tau_CaT_MSN(y(2));
       % mInf hInf
       mInf_CaT_MSN_Dendrite = infinite(y(2),CaT_MSN.mVhalf,CaT_MSN.mk);
       hInf_CaT_MSN_Dendrite = infinite(y(2),CaT_MSN.hVhalf,CaT_MSN.hk);
       
       %  49 50 m,h
       dy(49) = (mInf_CaT_MSN_Dendrite-y(49))/taum_CaT_MSN_Dendrite;
       dy(50) = (hInf_CaT_MSN_Dendrite-y(50))/tauh_CaT_MSN_Dendrite;
       
       g_CaT_MSN_Dendrite = g_Ca_MSN(y(2),y(51),Ca_0_MSN);
       I_CaT_MSN_Dendrite = g_CaT_MSN_Dendrite*CaT_MSN.Pbar*y(49)*y(49)*y(49)*y(50);%------mA
       
       %----------------------51 CaLdyn ----------------------------
       
       I_CaL_MSN_Dendrite = I_CaL1_2_MSN_Dendrite+I_CaL1_3_MSN_Dendrite+I_CaT_MSN_Dendrite;
       drive_channel_CaL_MSN_Dendrite = -drive_MSN * I_CaL_MSN_Dendrite /(2*F_MSN*depth_MSN);
       if drive_channel_CaL_MSN_Dendrite<=0
           drive_channel_CaL_MSN_Dendrite = 0;
       end
       drive_pump_CaL_MSN_Dendrite = -Kt_Ca_MSN * y(51)/(y(51)+Kd_Ca_MSN);
       dy(51) = drive_channel_CaL_MSN_Dendrite + pump_MSN*drive_pump_CaL_MSN_Dendrite +(CaL_i_inf_MSN-y(51))/tau_R_MSN;
       
       
       %----------------------52 Ca ----------------------------
       I_Ca_MSN_Dendrite = I_CaN_MSN_Dendrite+I_CaQ_MSN_Dendrite+I_CaR_MSN_Dendrite;
       drive_channel_Ca_MSN_Dendrite = -drive_MSN * I_Ca_MSN_Dendrite /(2*F_MSN*depth_MSN);
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
       
       
       
       %----------------------I_L----------------------------
       I_L_MSN_Dendrite = Leak_MSN.gmax*(y(2)-Leak_MSN.rev);
       

       %----------------------Synaptic current---------------------
       % NAc MSN:1 , PV: 99 CB:143
       % PFC Pyras:189+55*0~19; PV:1289+47*0~2; CB:1430+49*0~1
       % Inhibitory
       % 
       T_Presynapse_NAcPV = 1/(1+exp(-(y(99)-2)/5));%Transmitter Concentration
       dy(57) = GABAa_NAcInter2MSN.alpha*T_Presynapse_NAcPV*(1-y(57))-GABAa_NAcInter2MSN.beta*y(57);
       I_GABAa_NAcPV2MSN = y(57)*GABAa_NAcInter2MSN.gmax*(y(1)-GABAa_NAcInter2MSN.rev)*...
           Net_Strength(1,2);

       % CB 
       T_Presynapse_NAcCB = 1/(1+exp(-(y(143)-2)/5));%Transmitter Concentration
       dy(58) = GABAa_NAcInter2MSN.alpha*T_Presynapse_NAcCB*(1-y(58))-GABAa_NAcInter2MSN.beta*y(58);
       I_GABAa_NAcCB2MSN = y(58)*GABAa_NAcInter2MSN.gmax*(y(1)-GABAa_NAcInter2MSN.rev)*...
           Net_Strength(1,3);

       I_MSN_Inhibitory = I_GABAa_NAcPV2MSN+I_GABAa_NAcCB2MSN;
       % Excitatory
       % PFC Pyras:189+55*0~19; PV:1289+47*0~2; CB:1430+49*0~1
       NMDA_Mg_Block_MSN = 1/(1+exp(-(y(1)+15)/16.13));
       T_Presynapse_PFCPyra2MSN = zeros(20,1);
       I_AMPA_PFCPyra2MSN = zeros(20,1);
       I_NMDA_PFCPyra2MSN = zeros(20,1);
       
       
       for PFC_Pyra_Num = 0:19
           T_Presynapse_PFCPyra2MSN(PFC_Pyra_Num+1) = 1/(1+exp(-(y(189+55*PFC_Pyra_Num)/5)));
           dy(59+2*PFC_Pyra_Num) = AMPA_PFCPyra2MSN.alpha*T_Presynapse_PFCPyra2MSN(PFC_Pyra_Num+1)*...
               (1-y(59+2*PFC_Pyra_Num))-AMPA_PFCPyra2MSN.beta*y(59+2*PFC_Pyra_Num);%AMPA Gates
           I_AMPA_PFCPyra2MSN(PFC_Pyra_Num+1) = AMPA_PFCPyra2MSN.gmax*y(59+2*PFC_Pyra_Num)*...
               (y(2)-AMPA_PFCPyra2MSN.rev)*Net_Strength(1,4+PFC_Pyra_Num);%Postsynapse:Soma of Pyramidal
           dy(60+2*PFC_Pyra_Num) = NMDA_PFCPyra2MSN.alpha*T_Presynapse_PFCPyra2MSN(PFC_Pyra_Num+1)*...
               (1-y(60+2*PFC_Pyra_Num))-NMDA_PFCPyra2MSN.beta*y(60+2*PFC_Pyra_Num);%NMDA Gates
           I_NMDA_PFCPyra2MSN(PFC_Pyra_Num+1) = NMDA_PFCPyra2MSN.gmax*y(60+2*PFC_Pyra_Num)*...
               NMDA_Mg_Block_MSN*(y(2)-NMDA_PFCPyra2MSN.rev)*Net_Strength(1,4+PFC_Pyra_Num);
       end
       I_MSN_Excitatory = sum(I_AMPA_PFCPyra2MSN) + sum(I_NMDA_PFCPyra2MSN);
       %------------main function-------------
       I_ion_MSN_Soma = (I_NaF_MSN_Soma+I_NaP_MSN_Soma+I_KAf_MSN_Soma+I_KAs_MSN_Soma+I_KIR_MSN_Soma+...
           I_KRP_MSN_Soma+I_CaL1_2_MSN_Soma+I_CaL1_3_MSN_Soma+I_CaN_MSN_Soma+...
           I_CaQ_MSN_Soma+I_CaR_MSN_Soma+I_CaT_MSN_Soma+I_L_MSN_Soma+I_SKKCa_MSN_Soma+I_BKKCa_MSN_Soma)*1000;%1000
       I_ion_MSN_Dendrite = (I_NaF_MSN_Dendrite+I_NaP_MSN_Dendrite+I_KAf_MSN_Dendrite+I_KAs_MSN_Dendrite+I_KIR_MSN_Dendrite+...
           I_CaL1_2_MSN_Dendrite+I_CaL1_3_MSN_Dendrite+I_CaN_MSN_Dendrite+...
           I_CaQ_MSN_Dendrite+I_CaR_MSN_Dendrite+I_CaT_MSN_Dendrite+I_L_MSN_Dendrite+I_SKKCa_MSN_Dendrite+I_BKKCa_MSN_Dendrite)*1000;
       
       dy(1) = -I_ion_MSN_Soma - 10*g_Rall_MSN(1)*(y(1)-y(2)) + I_Ext_MSN_Soma_Real(floor(t/0.02)+1) - I_MSN_Inhibitory;
       dy(2) = -I_ion_MSN_Dendrite - g_Rall_MSN(2)*(y(2)-y(1)) + I_Ext_MSN_Dendrite_Real(floor(t/0.02)+1) - I_MSN_Excitatory;
       %% NAc, PV
       % PV:   99:Soma; 100-101: Ion Channels; 102:CB GABAa Soma;
       %       103-142: 20*Glu(AMPA+NMDA), Soma, PFC;
       [infm_NAc_Inter_PV_Na,tauh_NAc_Inter_PV_Na, infh_NAc_Inter_PV_Na] = inftau_PFC_Inter_PV_mh_Na(y(99));
       dy(100) = Na_PFC_Inter_PV.phi * (infh_NAc_Inter_PV_Na - y(100)) / tauh_NAc_Inter_PV_Na;
       I_NAc_Inter_PV_Na = Na_PFC_Inter_PV.gmax * (infm_NAc_Inter_PV_Na ^ 3) * y(100) * (y(99) - Na_PFC_Inter_PV.rev);

       [taun_NAc_Inter_PV_K, infn_NAc_Inter_PV_K] = inftau_PFC_Inter_PV_n_K(y(99));
       dy(101) = K_PFC_Inter_PV.phi * (infn_NAc_Inter_PV_K - y(101)) / taun_NAc_Inter_PV_K;
       I_NAc_Inter_PV_K = K_PFC_Inter_PV.gmax * (y(101) ^ 4) * (y(99) - K_PFC_Inter_PV.rev);
       
       I_NAc_Inter_PV_Leak = Leak_PFC_Inter_PV.gmax * (y(99) - Leak_PFC_Inter_PV.rev);

       %Inhibitory
       %CB
       %GABAa
       T_Presynapse_NAcCB2NAcPV = 1/(1+exp(-(y(143)-2)/5));%Transmitter Concentration
       dy(102) = GABAa_NAcInter2NAcPV.alpha*T_Presynapse_NAcCB2NAcPV*...
           (1-y(102))-GABAa_NAcInter2NAcPV.beta*y(102);
       I_GABAa_NAcCB2NAcPV = y(102)*...
           GABAa_NAcInter2NAcPV.gmax*(y(99)-GABAa_NAcInter2NAcPV.rev)*...
           Net_Strength(2,3);
       I_NAc_Inter_PV_Inhibitory = I_GABAa_NAcCB2NAcPV;
       %Excitatory
       NMDA_Mg_Block_NAc_Inter_PV = 1/(1+exp(-(y(99)+15)/16.13));
       T_Presynapse_PFCPyra2NAcPV = zeros(20,1);
       I_AMPA_PFCPyra2NAcPV = zeros(20,1);
       I_NMDA_PFCPyra2NAcPV = zeros(20,1);
       for PFC_Pyra_Num = 0:19
           T_Presynapse_PFCPyra2NAcPV(PFC_Pyra_Num+1) = 1/(1+exp(-(y(189+55*PFC_Pyra_Num)/5)));
           dy(103+2*PFC_Pyra_Num) = AMPA_PFCPyra2NAcPV.alpha*T_Presynapse_PFCPyra2NAcPV(PFC_Pyra_Num+1)*...
               (1-y(103+2*PFC_Pyra_Num))-AMPA_PFCPyra2NAcPV.beta*y(103+2*PFC_Pyra_Num);%AMPA Gates
           I_AMPA_PFCPyra2NAcPV(PFC_Pyra_Num+1) = AMPA_PFCPyra2NAcPV.gmax*y(103+2*PFC_Pyra_Num)*...
               (y(99)-AMPA_PFCPyra2NAcPV.rev)*Net_Strength(2,4+PFC_Pyra_Num);%Postsynapse:Soma of Pyramidal
           if I_AMPA_PFCPyra2NAcPV(PFC_Pyra_Num+1)>0
               I_AMPA_PFCPyra2NAcPV(PFC_Pyra_Num+1) = 0;
           end
           dy(104+2*PFC_Pyra_Num) = NMDA_PFCPyra2NAcPV.alpha*T_Presynapse_PFCPyra2NAcPV(PFC_Pyra_Num+1)*...
               (1-y(104+2*PFC_Pyra_Num))-NMDA_PFCPyra2NAcPV.beta*y(104+2*PFC_Pyra_Num);%NMDA Gates
           I_NMDA_PFCPyra2NAcPV(PFC_Pyra_Num+1) = NMDA_PFCPyra2NAcPV.gmax*y(104+2*PFC_Pyra_Num)*...
               NMDA_Mg_Block_NAc_Inter_PV*(y(99)-NMDA_PFCPyra2NAcPV.rev)*Net_Strength(2,4+PFC_Pyra_Num);
           if I_NMDA_PFCPyra2NAcPV(PFC_Pyra_Num+1)>0
               I_NMDA_PFCPyra2NAcPV(PFC_Pyra_Num+1) = 0;
           end
       end
       I_NAc_Inter_PV_Excitatory = sum(I_AMPA_PFCPyra2NAcPV) + sum(I_NMDA_PFCPyra2NAcPV);
       I_NAc_Inter_PV_Synapse = I_NAc_Inter_PV_Inhibitory + I_NAc_Inter_PV_Excitatory;
       
       %Total
       I_NAc_Inter_PV_Ion = I_NAc_Inter_PV_Na+I_NAc_Inter_PV_K+I_NAc_Inter_PV_Leak;
       dy(99) = - I_NAc_Inter_PV_Ion + I_Ext_NAc_Inter_PV_Real(floor(t/0.02)+1) - I_NAc_Inter_PV_Synapse;    
       %% NAc, CB
       % CB:   143:Soma, 144-147: Ion Channels; 148:PV GABAa Soma;
       %       149-188: 20*Glu(AMPA+NMDA), Soma, PFC;
       
       [infm_NAc_Inter_CB_Na, tauh_NAc_Inter_CB_Na, infh_NAc_Inter_CB_Na] = inftau_PFC_Inter_CB_mh_Na(y(143));
       dy(144) = Na_PFC_Inter_CB.phi * (infh_NAc_Inter_CB_Na - y(144)) / tauh_NAc_Inter_CB_Na;
       I_NAc_Inter_CB_Na = Na_PFC_Inter_CB.gmax * (infm_NAc_Inter_CB_Na ^ 3) * y(144) * (y(143) - Na_PFC_Inter_CB.rev);
       
       %'''----------------K 3 n--------------------'''
       
       [taun_NAc_Inter_CB_K, infn_NAc_Inter_CB_K] = inftau_PFC_Inter_CB_n_K(y(143));
       dy(145) = K_PFC_Inter_CB.phi * (infn_NAc_Inter_CB_K - y(145)) / taun_NAc_Inter_CB_K;
       I_NAc_Inter_CB_K = K_PFC_Inter_CB.gmax * (y(145) ^ 4) * (y(143) - K_PFC_Inter_CB.rev);
       
       %'''----------------Ca m--------------------'''
       infm_NAc_Inter_CB_Ca = 1 / (1 + exp(-(y(143) + 20) / 9));
       I_NAc_Inter_CB_Ca = Ca_PFC_Inter_CB.gmax * (infm_NAc_Inter_CB_Ca ^ 2) * (y(143) - Ca_PFC_Inter_CB.rev);
       
       
       %'''-----------------Ca dynamics 4--------------------'''
       alpha_Cadyn_Inter_CB = 0.002;
       tau_Cadyn_Inter_CB = 80;
       I_NAc_Inter_CB_Cadyn = I_NAc_Inter_CB_Ca;
       dy(146) = -alpha_Cadyn_Inter_CB * I_NAc_Inter_CB_Cadyn - y(146) / tau_Cadyn_Inter_CB;
       
       %'''-----------------KCa--------------------'''
       KD_Inter_CB = 30;
       I_NAc_Inter_CB_KCa = KCa_PFC_Inter_CB.gmax * (y(146) / (y(146) + KD_Inter_CB)) * (y(143) - KCa_PFC_Inter_CB.rev);
       
       %'''-----------------Ih H 5--------------------'''
       
       [tauh_NAc_Inter_CB_h, infh_NAc_Inter_CB_h] = inftau_PFC_Inter_CB_h_h(y(143));
       dy(147) = (infh_NAc_Inter_CB_h - y(147)) / tauh_NAc_Inter_CB_h;
       I_NAc_Inter_CB_h = h_PFC_Inter_CB.gmax * y(147) * (y(143) - h_PFC_Inter_CB.rev);
       
       %'''-----------------Cl--------------------'''
       
       I_NAc_Inter_CB_Leak = Leak_PFC_Inter_CB.gmax * (y(143) - Leak_PFC_Inter_CB.rev);
       
       % Inhibitory
       % PV
       % GABAa
       T_Presynapse_NAcPV2NAcCB = 1/(1+exp(-(y(99)-2)/5));%Transmitter Concentration
       dy(148) = GABAa_NAcInter2NAcCB.alpha*T_Presynapse_NAcPV2NAcCB*...
           (1-y(148))-GABAa_NAcInter2NAcCB.beta*y(148);
       I_GABAa_NAcPV2NAcCB = y(148)*...
           GABAa_NAcInter2NAcCB.gmax*(y(143)-GABAa_NAcInter2NAcCB.rev)*...
           Net_Strength(3,2);
       I_NAc_Inter_CB_Inhibitory = I_GABAa_NAcPV2NAcCB;
       % Excitatory
       NMDA_Mg_Block_NAc_Inter_CB = 1/(1+exp(-(y(143)+15)/16.13));
       T_Presynapse_PFCPyra2NAcCB = zeros(20,1);
       I_AMPA_PFCPyra2NAcCB = zeros(20,1);
       I_NMDA_PFCPyra2NAcCB = zeros(20,1);
       for PFC_Pyra_Num = 0:19
           T_Presynapse_PFCPyra2NAcCB(PFC_Pyra_Num+1) = 1/(1+exp(-(y(189+55*PFC_Pyra_Num)/5)));
           dy(149+2*PFC_Pyra_Num) = AMPA_PFCPyra2NAcCB.alpha*T_Presynapse_PFCPyra2NAcCB(PFC_Pyra_Num+1)*...
               (1-y(149+2*PFC_Pyra_Num))-AMPA_PFCPyra2NAcCB.beta*y(149+2*PFC_Pyra_Num);%AMPA Gates
           I_AMPA_PFCPyra2NAcCB(PFC_Pyra_Num+1) = AMPA_PFCPyra2NAcCB.gmax*y(149+2*PFC_Pyra_Num)*...
               (y(143)-AMPA_PFCPyra2NAcCB.rev)*Net_Strength(3,4+PFC_Pyra_Num);%Postsynapse:Soma of Pyramidal
           dy(150+2*PFC_Pyra_Num) = NMDA_PFCPyra2NAcCB.alpha*T_Presynapse_PFCPyra2NAcCB(PFC_Pyra_Num+1)*...
               (1-y(150+2*PFC_Pyra_Num))-NMDA_PFCPyra2NAcCB.beta*y(150+2*PFC_Pyra_Num);%NMDA Gates
           I_NMDA_PFCPyra2NAcCB(PFC_Pyra_Num+1) = NMDA_PFCPyra2NAcCB.gmax*y(150+2*PFC_Pyra_Num)*...
               NMDA_Mg_Block_NAc_Inter_CB*(y(143)-NMDA_PFCPyra2NAcCB.rev)*Net_Strength(3,4+PFC_Pyra_Num);
       end
       I_NAc_Inter_CB_Excitatory = sum(I_AMPA_PFCPyra2NAcCB) + sum(I_NMDA_PFCPyra2NAcCB);
       I_NAc_Inter_CB_Synapse = I_NAc_Inter_CB_Inhibitory + I_NAc_Inter_CB_Excitatory;
       
       
       %'''Total'''
       I_NAc_Inter_CB_Ion = I_NAc_Inter_CB_Na + I_NAc_Inter_CB_K + I_NAc_Inter_CB_Ca + ...
           I_NAc_Inter_CB_KCa + I_NAc_Inter_CB_h + I_NAc_Inter_CB_Leak;
       dy(143) = -I_NAc_Inter_CB_Ion + I_Ext_NAc_Inter_CB_Real(floor(t/0.02)+1) -I_NAc_Inter_CB_Synapse;
       %% PFC
       %% Pyra
       %Pyra:  189:Soma, 190:Proximal, 191:Distal; 192-200:Ion Channels;
       %       201-203 3*PV,GABAa Soma; 204-205 2*CB,GABAa
       %       206-243: 19*Pyras AMPA+NMDA
       %       In total: 20*55=1100, 189-1288: Pyras

       for k_PFC_Pyra = 0:19
           [infm_PFC_Pyra_Soma_Na, tauh_PFC_Pyra_Soma_Na, infh_PFC_Pyra_Soma_Na] = inftau_PFC_Pyra_Soma_mh_Na(y(189+55*k_PFC_Pyra));     
           dy(192+55*k_PFC_Pyra) = Na_PFC_Pyra_Soma.phi*(infh_PFC_Pyra_Soma_Na-y(192+55*k_PFC_Pyra))/tauh_PFC_Pyra_Soma_Na;
           I_PFC_Pyra_Soma_Na = Na_PFC_Pyra_Soma.gmax*(infm_PFC_Pyra_Soma_Na^3)*y(192+55*k_PFC_Pyra)*(y(189+55*k_PFC_Pyra)-Na_PFC_Pyra_Soma.rev);
           
           %-----------------Soma:K 5 n--------------------
           
           [taun_PFC_Pyra_Soma_K, infn_PFC_Pyra_Soma_K] = inftau_PFC_Pyra_Soma_n_K(y(189+55*k_PFC_Pyra));
           dy(193+55*k_PFC_Pyra) = K_PFC_Pyra_Soma.phi*(infn_PFC_Pyra_Soma_K-y(193+55*k_PFC_Pyra))/taun_PFC_Pyra_Soma_K;
           I_PFC_Pyra_Soma_K = K_PFC_Pyra_Soma.gmax*(y(193+55*k_PFC_Pyra)^4)*(y(189+55*k_PFC_Pyra)-K_PFC_Pyra_Soma.rev);
           
           %-----------------Soma:Ca m--------------------
           
           infm_PFC_Pyra_Soma_Ca = 1/(1+exp(-(y(189+55*k_PFC_Pyra)+20)/9));
           I_PFC_Pyra_Soma_Ca = Ca_PFC_Pyra_Soma.gmax(floor(t/0.02)+1)*(infm_PFC_Pyra_Soma_Ca^2)*(y(189+55*k_PFC_Pyra)-Ca_PFC_Pyra_Soma.rev);
           
           %-----------------Soma:CaN m 6--------------------
           
           alpham_PFC_Pyra_Soma_CaN = 0.0056;
           betam_PFC_Pyra_Soma_CaN = 0.002;
           infm_PFC_Pyra_Soma_CaN = alpham_PFC_Pyra_Soma_CaN*(y(195+55*k_PFC_Pyra)^2)/(alpham_PFC_Pyra_Soma_CaN*(y(195+55*k_PFC_Pyra)^2)+betam_PFC_Pyra_Soma_CaN);
           taum_PFC_Pyra_Soma_CaN = 1/(alpham_PFC_Pyra_Soma_CaN*(y(195+55*k_PFC_Pyra)^2)+betam_PFC_Pyra_Soma_CaN);
           dy(194+55*k_PFC_Pyra) = (infm_PFC_Pyra_Soma_CaN-y(194+55*k_PFC_Pyra))/taum_PFC_Pyra_Soma_CaN;
           I_PFC_Pyra_Soma_CaN = CaN_PFC_Pyra_Soma.gmax(floor(t/0.02)+1)*(y(194+55*k_PFC_Pyra)^2)*(y(189+55*k_PFC_Pyra)-CaN_PFC_Pyra_Soma.rev);
           
           %-----------------Soma:Ca dynamics 7--------------------
           
           alpha_PFC_Pyra_Soma_Cadyn = 0.000667;
           tau_PFC_Pyra_Soma_Cadyn = 240;
           I_PFC_Pyra_Soma_Cadyn = I_PFC_Pyra_Soma_Ca+I_PFC_Pyra_Soma_CaN;
           dy(195+55*k_PFC_Pyra) = -alpha_PFC_Pyra_Soma_Cadyn*I_PFC_Pyra_Soma_Cadyn-y(195+55*k_PFC_Pyra)/tau_PFC_Pyra_Soma_Cadyn;
           
           %-----------------Soma:Cl--------------------
           I_PFC_Pyra_Soma_Leak = Leak_PFC_Pyra.gmax*(y(189+55*k_PFC_Pyra)-Leak_PFC_Pyra.rev);
           
           
           
           %-----------Proximal Dendrites----------------
           %-----------------Proximal:NaP 8 h--------------------
           
           [infm_PFC_Pyra_Proximal_Na,tauh_PFC_Pyra_Proximal_Na, infh_PFC_Pyra_Proximal_Na] = inftau_PFC_Pyra_Proximal_mh_Na(y(190+55*k_PFC_Pyra),Dopamine_Ratio(floor(t/0.02)+1));
           dy(196+55*k_PFC_Pyra) = (infh_PFC_Pyra_Proximal_Na-y(196+55*k_PFC_Pyra))/tauh_PFC_Pyra_Proximal_Na;
           I_PFC_Pyra_Proximal_NaP = NaP_PFC_Pyra_Proximal.gmax*(infm_PFC_Pyra_Proximal_Na^3)*y(196+55*k_PFC_Pyra)*(y(190+55*k_PFC_Pyra)-NaP_PFC_Pyra_Proximal.rev);
           
           %-----------------Proximal:KS q r 9 10--------------------
           
           %        16;
           [tauq_PFC_Pyra_Proximal_KS, infq_PFC_Pyra_Proximal_KS,taur_PFC_Pyra_Proximal_KS, infr_PFC_Pyra_Proximal_KS] = inftau_PFC_Pyra_Proximal_qr_KS(y(190+55*k_PFC_Pyra));
           dy(197+55*k_PFC_Pyra) = (infq_PFC_Pyra_Proximal_KS-y(197+55*k_PFC_Pyra))/tauq_PFC_Pyra_Proximal_KS;
           dy(198+55*k_PFC_Pyra) = (infr_PFC_Pyra_Proximal_KS-y(198+55*k_PFC_Pyra))/taur_PFC_Pyra_Proximal_KS;
           I_PFC_Pyra_Proximal_KS = KS_PFC_Pyra_Proximal.gmax(floor(t/0.02)+1)*y(197+55*k_PFC_Pyra)*y(198+55*k_PFC_Pyra)*(y(190+55*k_PFC_Pyra)-KS_PFC_Pyra_Proximal.rev);
           
           
           %-----------------Proximal:Cl--------------------
           
           I_PFC_Pyra_Proximal_Leak =  Leak_PFC_Pyra.gmax*(y(190+55*k_PFC_Pyra)-Leak_PFC_Pyra.rev);
           
           
           %-----------Distal Dendrites----------------
           %-----------------Distal:A a b 11 12--------------------
           
           
           [taua_PFC_Pyra_Distal_A, infa_PFC_Pyra_Distal_A,taub_PFC_Pyra_Distal_A, infb_PFC_Pyra_Distal_A] = inftau_PFC_Pyra_Distal_ab_A(y(191+55*k_PFC_Pyra));
           dy(199+55*k_PFC_Pyra) = (infa_PFC_Pyra_Distal_A-y(199+55*k_PFC_Pyra))/taua_PFC_Pyra_Distal_A;
           dy(200+55*k_PFC_Pyra) = (infb_PFC_Pyra_Distal_A-y(200+55*k_PFC_Pyra))/taub_PFC_Pyra_Distal_A;
           I_PFC_Pyra_Distal_A = A_PFC_Pyra_Distal.gmax*(y(199+55*k_PFC_Pyra)^4)*y(200+55*k_PFC_Pyra)*(y(191+55*k_PFC_Pyra)-A_PFC_Pyra_Distal.rev);
           
           %-----------------Distal:Ca--------------------
           infm_PFC_Pyra_Distal_Ca = 1/(1+exp(-(y(191+55*k_PFC_Pyra)+20)/9));
           I_PFC_Pyra_Distal_Ca = Ca_PFC_Pyra_Distal.gmax(floor(t/0.02)+1)*(infm_PFC_Pyra_Distal_Ca^2)*(y(191+55*k_PFC_Pyra)-Ca_PFC_Pyra_Distal.rev);
           
           
           %-----------------Distal:Cl--------------------
           
           I_PFC_Pyra_Distal_Leak =  Leak_PFC_Pyra.gmax*(y(191+55*k_PFC_Pyra)-Leak_PFC_Pyra.rev);
           T_Presynapse_PFCPV2PFCPyra = zeros(3,1);
           I_GABAa_PFCPV2PFCPyra = zeros(3,1);
           % Inhibitory
           for k_PFC_PV_Pre = 0:2
               T_Presynapse_PFCPV2PFCPyra(k_PFC_PV_Pre+1) = 1/(1+exp(-(y(1289+47*k_PFC_PV_Pre)-2)/5));%Transmitter Concentration
               %GABAa
               dy(201+55*k_PFC_Pyra+k_PFC_PV_Pre) = GABAa_PFCInter2PFCPyra.alpha*T_Presynapse_PFCPV2PFCPyra(k_PFC_PV_Pre+1)*...
                   (1-y(201+55*k_PFC_Pyra+k_PFC_PV_Pre))-GABAa_PFCInter2PFCPyra.beta*y(201+55*k_PFC_Pyra+k_PFC_PV_Pre);
               I_GABAa_PFCPV2PFCPyra(k_PFC_PV_Pre+1) = y(201+55*k_PFC_Pyra+k_PFC_PV_Pre)*...
                   GABAa_PFCInter2PFCPyra.gmax*(y(189+55*k_PFC_Pyra)-GABAa_PFCInter2PFCPyra.rev)*...
                   Net_Strength(4+k_PFC_Pyra,24+k_PFC_PV_Pre);
           end
           % 626~631CB Soma: 2070 + 57*(0~1);
           T_Presynapse_PFCCB2PFCPyra = zeros(2,1);
           I_GABAa_PFCCB2PFCPyra = zeros(2,1);
           for k_PFC_CB_Pre = 0:1
               T_Presynapse_PFCCB2PFCPyra(k_PFC_CB_Pre+1) = 1/(1+exp(-(y(1430+49*k_PFC_CB_Pre)-2)/5));%Transmitter Concentration
               %GABAa
               dy(204+55*k_PFC_Pyra+k_PFC_CB_Pre) = GABAa_PFCInter2PFCPyra.alpha*T_Presynapse_PFCCB2PFCPyra(k_PFC_CB_Pre+1)*...
                   (1-y(204+55*k_PFC_Pyra+k_PFC_CB_Pre))-GABAa_PFCInter2PFCPyra.beta*y(204+55*k_PFC_Pyra+k_PFC_CB_Pre);
               I_GABAa_PFCCB2PFCPyra(k_PFC_CB_Pre+1) = y(204+55*k_PFC_Pyra+k_PFC_CB_Pre)*...
                   GABAa_PFCInter2PFCPyra.gmax*(y(189+55*k_PFC_Pyra)-GABAa_PFCInter2PFCPyra.rev)*...
                   Net_Strength(4+k_PFC_Pyra,27+k_PFC_CB_Pre);
           end
           I_PFC_Pyra_Inhibitory = sum(I_GABAa_PFCPV2PFCPyra)+sum(I_GABAa_PFCCB2PFCPyra);
           % Pyras
           NMDA_Mg_Block_PFC_Pyra = 1/(1+exp(-(y(189+55*k_PFC_Pyra)+15)/16.13));
           k_PFC_Pyra_Pre_No = 0; %count
           T_Presynapse_PFCPyra2PFCPyra = zeros(19,1);
           I_AMPA_PFCPyra2PFCPyra = zeros(19,1);
           I_NMDA_PFCPyra2PFCPyra = zeros(19,1);
           for k_PFC_Pyra_Pre = 0:19
               if k_PFC_Pyra_Pre == k_PFC_Pyra
                   continue
               else
                   T_Presynapse_PFCPyra2PFCPyra(k_PFC_Pyra_Pre_No+1) = 1/(1+exp(-(y(189+55*k_PFC_Pyra_Pre)/5)));
                   dy(206+k_PFC_Pyra_Pre_No*2+55*k_PFC_Pyra) = AMPA_PFCPyra2PFCPyra.alpha*T_Presynapse_PFCPyra2PFCPyra(k_PFC_Pyra_Pre_No+1)*...
                       (1-y(206+k_PFC_Pyra_Pre_No*2+55*k_PFC_Pyra))-AMPA_PFCPyra2PFCPyra.beta*y(206+k_PFC_Pyra_Pre_No*2+55*k_PFC_Pyra);%AMPA Gates
                   I_AMPA_PFCPyra2PFCPyra(k_PFC_Pyra_Pre_No+1) = AMPA_PFCPyra2PFCPyra.gmax*y(206+k_PFC_Pyra_Pre_No*2+55*k_PFC_Pyra)*...
                       (y(191 + 55*k_PFC_Pyra)-AMPA_PFCPyra2PFCPyra.rev)*Net_Strength(4+k_PFC_Pyra,4+k_PFC_Pyra_Pre);%Postsynapse:Soma of Pyramidal
                   dy(207+k_PFC_Pyra_Pre_No*2+55*k_PFC_Pyra) = NMDA_PFCPyra2PFCPyra.alpha*T_Presynapse_PFCPyra2PFCPyra(k_PFC_Pyra_Pre_No+1)*...
                       (1-y(207+k_PFC_Pyra_Pre_No*2+55*k_PFC_Pyra))-NMDA_PFCPyra2PFCPyra.beta*y(207+k_PFC_Pyra_Pre_No*2+55*k_PFC_Pyra);%NMDA Gates
                   I_NMDA_PFCPyra2PFCPyra(k_PFC_Pyra_Pre_No+1) = NMDA_PFCPyra2PFCPyra.gmax*y(207+k_PFC_Pyra_Pre_No*2+55*k_PFC_Pyra)*...
                       NMDA_Mg_Block_PFC_Pyra*(y(191 + 55*k_PFC_Pyra)-NMDA_PFCPyra2PFCPyra.rev)*Net_Strength(4+k_PFC_Pyra,4+k_PFC_Pyra_Pre);
                   k_PFC_Pyra_Pre_No = k_PFC_Pyra_Pre_No+1;
               end
           end
           I_PFC_Pyra_Excitatory = sum(I_AMPA_PFCPyra2PFCPyra)+sum(I_NMDA_PFCPyra2PFCPyra); 
           %-----------------Total--------------------
           
           I_PFC_Pyra_Soma_Ion = I_PFC_Pyra_Soma_Na+I_PFC_Pyra_Soma_K+I_PFC_Pyra_Soma_Ca+I_PFC_Pyra_Soma_CaN+I_PFC_Pyra_Soma_Leak;
           %         I_Soma = I_Na_Soma+I_K_Soma;
           gc1 = 0.75;
           p1 = 0.5; 
           dy(189 + 55*k_PFC_Pyra) = -I_PFC_Pyra_Soma_Ion-gc1*(y(189 + 55*k_PFC_Pyra)-y(190 + 55*k_PFC_Pyra))/p1+I_Ext_PFC_Pyra_Soma_Real(floor(t/0.02)+1)-I_PFC_Pyra_Inhibitory;
           I_PFC_Pyra_Proximal_Ion = I_PFC_Pyra_Proximal_NaP+I_PFC_Pyra_Proximal_KS+I_PFC_Pyra_Proximal_Leak;
           gc2 = 0.25;
           p2 = 0.3;
           dy(190 + 55*k_PFC_Pyra) = -I_PFC_Pyra_Proximal_Ion-gc1*(y(190 + 55*k_PFC_Pyra)-y(189 + 55*k_PFC_Pyra))/p2 -gc2*(y(190 + 55*k_PFC_Pyra)-y(191 + 55*k_PFC_Pyra))/p2+I_Ext_PFC_Pyra_Proximal_Real(floor(t/0.02)+1);
           I_PFC_Pyra_Distal_Ion = I_PFC_Pyra_Distal_A+I_PFC_Pyra_Distal_Ca+I_PFC_Pyra_Distal_Leak;
           %         I_Distal = I_A_Distal+I_L_Distal;
           dy(191 + 55*k_PFC_Pyra) = -I_PFC_Pyra_Distal_Ion-gc2*(y(191 + 55*k_PFC_Pyra)-y(190 + 55*k_PFC_Pyra))/(1-p1-p2)+I_Ext_PFC_Pyra_Distal_Real(floor(t/0.02)+1)-I_PFC_Pyra_Excitatory;
       end
       %% PV
       %PV:    1289:Soma; 1290-1291: Ion Channels; 1292-1293: 2*PV, GABAa Soma;
       %       1294-1295: 2*CB, GABAa Soma; 1296-1335: 20*PV,Glus;
       %       In total: 3*47 = 141, 1289-1429,PVs
       for k_PFC_PV = 0:2
           [infm_PFC_Inter_PV_Na,tauh_PFC_Inter_PV_Na, infh_PFC_Inter_PV_Na] = inftau_PFC_Inter_PV_mh_Na(y(1289+47*k_PFC_PV));
           dy(1290+47*k_PFC_PV) = Na_PFC_Inter_PV.phi * (infh_PFC_Inter_PV_Na - y(1290+47*k_PFC_PV)) / tauh_PFC_Inter_PV_Na;
           I_PFC_Inter_PV_Na = Na_PFC_Inter_PV.gmax * (infm_PFC_Inter_PV_Na ^ 3) * y(1290+47*k_PFC_PV) * (y(1289+47*k_PFC_PV) - Na_PFC_Inter_PV.rev);
           
           [taun_PFC_Inter_PV_K, infn_PFC_Inter_PV_K] = inftau_PFC_Inter_PV_n_K(y(1289+47*k_PFC_PV));
           dy(1291+47*k_PFC_PV) = K_PFC_Inter_PV.phi * (infn_PFC_Inter_PV_K - y(1291+47*k_PFC_PV)) / taun_PFC_Inter_PV_K;
           I_PFC_Inter_PV_K = K_PFC_Inter_PV.gmax * (y(1291+47*k_PFC_PV) ^ 4) * (y(1289+47*k_PFC_PV) - K_PFC_Inter_PV.rev);
           
           I_PFC_Inter_PV_Leak = Leak_PFC_Inter_PV.gmax * (y(1289+47*k_PFC_PV) - Leak_PFC_Inter_PV.rev);
           
           %Inhibitory
           %PV
           T_Presynapse_PFCPV2PFCPV = zeros(2,1);
           I_GABAa_PFCPV2PFCPV = zeros(2,1);
           k_PFC_PV_Pre_No = 0; %count
           for k_PFC_PV_Pre = 0:2
               if k_PFC_PV_Pre == k_PFC_PV
                   continue
               else
                   T_Presynapse_PFCPV2PFCPV(k_PFC_PV_Pre_No+1) = 1/(1+exp(-(y(1289+47*k_PFC_PV_Pre)/5)));
                   %GABAa
                   dy(1292+47*k_PFC_PV+k_PFC_PV_Pre_No) = GABAa_PFCInter2PFCPV.alpha*T_Presynapse_PFCPV2PFCPV(k_PFC_PV_Pre_No+1)*...
                       (1-y(1292+47*k_PFC_PV+k_PFC_PV_Pre_No))-GABAa_PFCInter2PFCPV.beta*y(1292+47*k_PFC_PV+k_PFC_PV_Pre_No);
                   I_GABAa_PFCPV2PFCPV(k_PFC_PV_Pre_No+1) = y(1292+47*k_PFC_PV+k_PFC_PV_Pre_No)*...
                       GABAa_PFCInter2PFCPV.gmax*(y(1289+47*k_PFC_PV)-GABAa_PFCInter2PFCPV.rev)*...
                       Net_Strength(24+k_PFC_PV,24+k_PFC_PV_Pre);
                   k_PFC_PV_Pre_No = k_PFC_PV_Pre_No+1; 
               end
           end  
           % CB
           T_Presynapse_PFCCB2PFCPV = zeros(2,1);
           I_GABAa_PFCCB2PFCPV = zeros(2,1);
           for k_PFC_CB_Pre = 0:1
               T_Presynapse_PFCCB2PFCPV(k_PFC_CB_Pre+1) = 1/(1+exp(-(y(1430+49*k_PFC_CB_Pre)-2)/5));%Transmitter Concentration
               %GABAa
               dy(1294+k_PFC_CB_Pre+47*k_PFC_PV) = GABAa_PFCInter2PFCPV.alpha*T_Presynapse_PFCCB2PFCPV(k_PFC_CB_Pre+1)*...
                   (1-y(1294+k_PFC_CB_Pre+47*k_PFC_PV))-GABAa_PFCInter2PFCPV.beta*y(1294+k_PFC_CB_Pre+47*k_PFC_PV);
               I_GABAa_PFCCB2PFCPV(k_PFC_CB_Pre+1) = y(1294+k_PFC_CB_Pre+47*k_PFC_PV)*...
                   GABAa_PFCInter2PFCPV.gmax*(y(1289+47*k_PFC_PV)-GABAa_PFCInter2PFCPV.rev)*...
                   Net_Strength(24+k_PFC_PV,27+k_PFC_CB_Pre);
           end
           I_PFC_Inter_PV_Inhibitory = sum(I_GABAa_PFCPV2PFCPV)+sum(I_GABAa_PFCCB2PFCPV);
           %Excitatory
           % 632~669 65*(0~19): 19*Pyras,AMPA+NMDA;
           NMDA_Mg_Block_PFC_PV = 1/(1+exp(-(y(1289+47*k_PFC_PV)+15)/16.13));
           T_Presynapse_PFCPyra2PFCPV = zeros(20,1);
           I_AMPA_PFCPyra2PFCPV = zeros(20,1);
           I_NMDA_PFCPyra2PFCPV = zeros(20,1);
           for k_PFC_Pyra_Pre = 0:19
               T_Presynapse_PFCPyra2PFCPV(k_PFC_Pyra_Pre+1) = 1/(1+exp(-(y(189+55*k_PFC_Pyra_Pre)/5)));
               dy(1296+47*k_PFC_PV+k_PFC_Pyra_Pre*2) = AMPA_PFCPyra2PFCPV.alpha*T_Presynapse_PFCPyra2PFCPV(k_PFC_Pyra_Pre+1)*...
                   (1-y(1296+47*k_PFC_PV+k_PFC_Pyra_Pre*2))-AMPA_PFCPyra2PFCPV.beta*y(1296+47*k_PFC_PV+k_PFC_Pyra_Pre*2);%AMPA Gates
               I_AMPA_PFCPyra2PFCPV(k_PFC_Pyra_Pre+1) = AMPA_PFCPyra2PFCPV.gmax*y(1296+47*k_PFC_PV+k_PFC_Pyra_Pre*2)*...
                   (y(1289+47*k_PFC_PV)-AMPA_PFCPyra2PFCPV.rev)*Net_Strength(24+k_PFC_PV,4+k_PFC_Pyra_Pre);%Postsynapse:Soma of Pyramidal
               dy(1297+47*k_PFC_PV+k_PFC_Pyra_Pre*2) = NMDA_PFCPyra2PFCPV.alpha*T_Presynapse_PFCPyra2PFCPV(k_PFC_Pyra_Pre+1)*...
                   (1-y(1297+47*k_PFC_PV+k_PFC_Pyra_Pre*2))-NMDA_PFCPyra2PFCPV.beta*y(1297+47*k_PFC_PV+k_PFC_Pyra_Pre*2);%NMDA Gates
               I_NMDA_PFCPyra2PFCPV(k_PFC_Pyra_Pre+1) = NMDA_PFCPyra2PFCPV.gmax*y(1297+47*k_PFC_PV+k_PFC_Pyra_Pre*2)*...
                   NMDA_Mg_Block_PFC_PV*(y(1289+47*k_PFC_PV)-NMDA_PFCPyra2PFCPV.rev)*Net_Strength(24+k_PFC_PV,4+k_PFC_Pyra_Pre);
           end
           I_PFC_Inter_PV_Excitatory = sum(I_AMPA_PFCPyra2PFCPV)+sum(I_NMDA_PFCPyra2PFCPV);
           I_PFC_Inter_PV_Synapse = I_PFC_Inter_PV_Inhibitory+I_PFC_Inter_PV_Excitatory;
           
           %Total
           I_PFC_Inter_PV_Ion = I_PFC_Inter_PV_Na+I_PFC_Inter_PV_K+I_PFC_Inter_PV_Leak;
           dy(1289+47*k_PFC_PV) = - I_PFC_Inter_PV_Ion + I_Ext_PFC_Inter_PV_Real(floor(t/0.02)+1) - I_PFC_Inter_PV_Synapse;
       end
       %% CB
       %CB:    1430:Soma, 1431-1434: Ion Channels, 1435-1437: 3*PV, GABAa Soma;
       %       1438: 1*CB, GABAa Soma; 1439-1478: 20*PV,Glus;
       %       In total: 2*49 = 98, 1430-1527: CBs
       for k_PFC_CB = 0:1
           [infm_PFC_Inter_CB_Na, tauh_PFC_Inter_CB_Na, infh_PFC_Inter_CB_Na] = inftau_PFC_Inter_CB_mh_Na(y(1430+49*k_PFC_CB));
           dy(1431+49*k_PFC_CB) = Na_PFC_Inter_CB.phi * (infh_PFC_Inter_CB_Na - y(1431+49*k_PFC_CB)) / tauh_PFC_Inter_CB_Na;
           I_PFC_Inter_CB_Na = Na_PFC_Inter_CB.gmax * (infm_PFC_Inter_CB_Na ^ 3) * y(1431+49*k_PFC_CB) * (y(1430+49*k_PFC_CB) - Na_PFC_Inter_CB.rev);
           
           %'''----------------K 3 n--------------------'''
           
           [taun_PFC_Inter_CB_K, infn_PFC_Inter_CB_K] = inftau_PFC_Inter_CB_n_K(y(1430+49*k_PFC_CB));
           dy(1432+49*k_PFC_CB) = K_PFC_Inter_CB.phi * (infn_PFC_Inter_CB_K - y(1432+49*k_PFC_CB)) / taun_PFC_Inter_CB_K;
           I_PFC_Inter_CB_K = K_PFC_Inter_CB.gmax * (y(1432+49*k_PFC_CB) ^ 4) * (y(1430+49*k_PFC_CB) - K_PFC_Inter_CB.rev);
           
           %'''----------------Ca m--------------------'''
           infm_PFC_Inter_CB_Ca = 1 / (1 + exp(-(y(1430+49*k_PFC_CB) + 20) / 9));
           I_PFC_Inter_CB_Ca = Ca_PFC_Inter_CB.gmax * (infm_PFC_Inter_CB_Ca ^ 2) * (y(1430+49*k_PFC_CB) - Ca_PFC_Inter_CB.rev);
           
           
           %'''-----------------Ca dynamics 4--------------------'''
           alpha_Cadyn_Inter_CB = 0.002;
           tau_Cadyn_Inter_CB = 80;
           I_PFC_Inter_CB_Cadyn = I_PFC_Inter_CB_Ca;
           dy(1433+49*k_PFC_CB) = -alpha_Cadyn_Inter_CB * I_PFC_Inter_CB_Cadyn - y(1433+49*k_PFC_CB) / tau_Cadyn_Inter_CB;
           
           %'''-----------------KCa--------------------'''
           KD_Inter_CB = 30;
           I_PFC_Inter_CB_KCa = KCa_PFC_Inter_CB.gmax * (y(1433+49*k_PFC_CB) / (y(1433+49*k_PFC_CB) + KD_Inter_CB)) * (y(1430+49*k_PFC_CB) - KCa_PFC_Inter_CB.rev);
           
           %'''-----------------Ih H 5--------------------'''
           
           [tauh_PFC_Inter_CB_h, infh_PFC_Inter_CB_h] = inftau_PFC_Inter_CB_h_h(y(1430+49*k_PFC_CB));
           dy(1434+49*k_PFC_CB) = (infh_PFC_Inter_CB_h - y(1434+49*k_PFC_CB)) / tauh_PFC_Inter_CB_h;
           I_PFC_Inter_CB_h = h_PFC_Inter_CB.gmax * y(1434+49*k_PFC_CB) * (y(1430+49*k_PFC_CB) - h_PFC_Inter_CB.rev);
           
           %'''-----------------Cl--------------------'''
           
           I_PFC_Inter_CB_Leak = Leak_PFC_Inter_CB.gmax * (y(1430+49*k_PFC_CB) - Leak_PFC_Inter_CB.rev);
           
           % Inhibitory
           % PV
           % GABAa
           T_Presynapse_PFCPV2PFCCB = zeros(3,1);
           I_GABAa_PFCPV2PFCCB = zeros(3,1);
           for k_PFC_PV_Pre = 0:2
               T_Presynapse_PFCPV2PFCCB(k_PFC_PV_Pre+1) = 1/(1+exp(-(y(1289+47*k_PFC_PV_Pre)-2)/5));%Transmitter Concentration
               %GABAa
               dy(1435+49*k_PFC_CB+k_PFC_PV_Pre) = GABAa_PFCInter2PFCCB.alpha*T_Presynapse_PFCPV2PFCCB(k_PFC_PV_Pre+1)*...
                   (1-y(1435+49*k_PFC_CB+k_PFC_PV_Pre))-GABAa_PFCInter2PFCCB.beta*y(1435+49*k_PFC_CB+k_PFC_PV_Pre);
               I_GABAa_PFCPV2PFCCB(k_PFC_PV_Pre+1) = y(1435+49*k_PFC_CB+k_PFC_PV_Pre)*...
                   GABAa_PFCInter2PFCCB.gmax*(y(1430+49*k_PFC_CB)-GABAa_PFCInter2PFCCB.rev)*...
                   Net_Strength(27+k_PFC_CB,24+k_PFC_PV_Pre);
           end
           
           %CB
           for k_PFC_CB_Pre = 0:1
               if k_PFC_CB_Pre == k_PFC_CB
                   continue
               else
                   T_Presynapse_PFCCB2PFCCB = 1/(1+exp(-(y(1430+49*k_PFC_CB_Pre)-2)/5));%Transmitter Concentration
                   %GABAa
                   dy(1438+49*k_PFC_CB) = GABAa_PFCInter2PFCCB.alpha*T_Presynapse_PFCCB2PFCCB*...
                       (1-y(1438+49*k_PFC_CB))-GABAa_PFCInter2PFCCB.beta*y(1438+49*k_PFC_CB);
                   I_GABAa_PFCCB2PFCCB = y(1438+49*k_PFC_CB)*...
                       GABAa_PFCInter2PFCCB.gmax*(y(1430+49*k_PFC_CB)-GABAa_PFCInter2PFCCB.rev)*...
                       Net_Strength(27+k_PFC_CB,27+k_PFC_CB_Pre);
               end
           end
           I_PFC_Inter_CB_Inhibitory = sum(I_GABAa_PFCPV2PFCCB)+I_GABAa_PFCCB2PFCCB;
           
           %Excitatory
           %Pyras
           T_Presynapse_PFCPyra2PFCCB = zeros(20,1);
           I_AMPA_PFCPyra2PFCCB = zeros(20,1);
           I_NMDA_PFCPyra2PFCCB = zeros(20,1);
           NMDA_Mg_Block_PFC_CB = 1/(1+exp(-(y(1430+49*k_PFC_CB)+15)/16.13));
           for k_PFC_Pyra_Pre = 0:19
               T_Presynapse_PFCPyra2PFCCB(k_PFC_Pyra_Pre+1) = 1/(1+exp(-(y(189+55*k_PFC_Pyra_Pre)/5)));
               dy(1439+49*k_PFC_CB+2*k_PFC_Pyra_Pre) = AMPA_PFCPyra2PFCCB.alpha*T_Presynapse_PFCPyra2PFCCB(k_PFC_Pyra_Pre+1)*...
                   (1-y(1439+49*k_PFC_CB+2*k_PFC_Pyra_Pre))-AMPA_PFCPyra2PFCCB.beta*y(1439+49*k_PFC_CB+2*k_PFC_Pyra_Pre);%AMPA Gates
               I_AMPA_PFCPyra2PFCCB(k_PFC_Pyra_Pre+1) = AMPA_PFCPyra2PFCCB.gmax*y(1439+49*k_PFC_CB+2*k_PFC_Pyra_Pre)*...
                   (y(1430+49*k_PFC_CB)-AMPA_PFCPyra2PFCCB.rev)*Net_Strength(27+k_PFC_CB,4+k_PFC_Pyra_Pre);%Postsynapse:Soma of Pyramidal
               dy(1440+49*k_PFC_CB+2*k_PFC_Pyra_Pre) = NMDA_PFCPyra2PFCCB.alpha*T_Presynapse_PFCPyra2PFCCB(k_PFC_Pyra_Pre+1)*...
                   (1-y(1440+49*k_PFC_CB+2*k_PFC_Pyra_Pre))-NMDA_PFCPyra2PFCCB.beta*y(1440+49*k_PFC_CB+2*k_PFC_Pyra_Pre);%NMDA Gates
               I_NMDA_PFCPyra2PFCCB(k_PFC_Pyra_Pre+1) = NMDA_PFCPyra2PFCCB.gmax*y(1440+49*k_PFC_CB+2*k_PFC_Pyra_Pre)*...
                   NMDA_Mg_Block_PFC_CB*(y(1430+49*k_PFC_CB)-NMDA_PFCPyra2PFCCB.rev)*Net_Strength(27+k_PFC_CB,4+k_PFC_Pyra_Pre);
           end
           I_PFC_Inter_CB_Excitatory = sum(I_AMPA_PFCPyra2PFCCB)+sum(I_NMDA_PFCPyra2PFCCB);
           I_PFC_Inter_CB_Synapse = I_PFC_Inter_CB_Inhibitory+I_PFC_Inter_CB_Excitatory;
           
           
           
           
           %'''Total'''
           I_PFC_Inter_CB_Ion = I_PFC_Inter_CB_Na + I_PFC_Inter_CB_K + I_PFC_Inter_CB_Ca + ...
               I_PFC_Inter_CB_KCa + I_PFC_Inter_CB_h + I_PFC_Inter_CB_Leak;
           dy(1430+49*k_PFC_CB) = -I_PFC_Inter_CB_Ion + I_Ext_PFC_Inter_CB_Real(floor(t/0.02)+1) -I_PFC_Inter_CB_Synapse;
           
       end
    end
end

