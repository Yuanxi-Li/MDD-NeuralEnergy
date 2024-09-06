function [Small_Range,Medium_Range,Large_Range,Full_Range,Small_Data,Medium_Data,Large_Data,Full_Data] = test_FiringSequence_Dopamineratio_Change_2
% Experiment 1
I_Ext_MSN_Soma=3;
I_Ext_MSN_Dendrite=3;
I_Ext_PFC_Pyra_Soma=2;
I_Ext_PFC_Inter_CB=0.1;
I_Ext_PFC_Inter_PV=0.5;
I_Ext_NAc_Inter_CB=0.1;
I_Ext_NAc_Inter_PV=0.5;

Seed_Ext = 26;


Small_Range = {};
Medium_Range = {};
Large_Range = {};
Full_Range = {};
Small_Data = {};
Medium_Data = {};
Large_Data = {};
Full_Data = {};

for Dopamine_Ratio_Range = 0:0.25:0.75
    t_NAC_MSN_Temp = {};
    t_NAC_PV_Temp = {};
    t_NAC_CB_Temp = {};
    t_PFC_Pyra_Temp = {};
    t_PFC_PV_Temp = {};
    t_PFC_CB_Temp = {};
    for Seed_Net = 1:3  
        for Seed_Dopamine = 51:1:53
            rng(Seed_Dopamine)
            Dopamine_Ratio = randi([Dopamine_Ratio_Range*1e6 (Dopamine_Ratio_Range+0.25)*1e6],1,3005/0.02)/1e6;            
            [t,y,t_NAC_MSN,t_NAC_PV,t_NAC_CB,t_PFC_Pyra,t_PFC_PV,t_PFC_CB] = FiringSequence_Dopamineratio_Change(Seed_Net, Seed_Ext, Dopamine_Ratio,...
                1*(I_Ext_MSN_Soma-1),1*(I_Ext_MSN_Dendrite-1),1*(I_Ext_PFC_Pyra_Soma-1),1*(I_Ext_PFC_Pyra_Soma-1),1*(I_Ext_PFC_Pyra_Soma-1),...
                I_Ext_PFC_Inter_CB,I_Ext_PFC_Inter_PV, I_Ext_NAc_Inter_CB, I_Ext_NAc_Inter_PV);
            t_NAC_MSN_Temp{Seed_Net,Seed_Dopamine-50} = t_NAC_MSN;
            t_NAC_PV_Temp{Seed_Net,Seed_Dopamine-50} = t_NAC_PV;
            t_NAC_CB_Temp{Seed_Net,Seed_Dopamine-50} = t_NAC_CB;
            t_PFC_Pyra_Temp{Seed_Net,Seed_Dopamine-50} = t_PFC_Pyra;
            t_PFC_PV_Temp{Seed_Net,Seed_Dopamine-50} = t_PFC_PV;
            t_PFC_CB_Temp{Seed_Net,Seed_Dopamine-50} = t_PFC_CB;
            y_Temp{Seed_Net,Seed_Dopamine-50} = y;
        end
    end
    if Dopamine_Ratio_Range == 0
        Small_Range = {t_NAC_MSN_Temp, t_NAC_PV_Temp, t_NAC_CB_Temp, t_PFC_Pyra_Temp, t_PFC_PV_Temp, t_PFC_CB_Temp};
        Small_Data = {y_Temp};
    elseif Dopamine_Ratio_Range == 0.25
        Medium_Range = {t_NAC_MSN_Temp, t_NAC_PV_Temp, t_NAC_CB_Temp, t_PFC_Pyra_Temp, t_PFC_PV_Temp, t_PFC_CB_Temp};
        Medium_Data = {y_Temp};
    elseif Dopamine_Ratio_Range == 0.5
        Large_Range = {t_NAC_MSN_Temp, t_NAC_PV_Temp, t_NAC_CB_Temp, t_PFC_Pyra_Temp, t_PFC_PV_Temp, t_PFC_CB_Temp};
        Large_Data = {y_Temp};
    elseif Dopamine_Ratio_Range == 0.75
        Full_Range = {t_NAC_MSN_Temp, t_NAC_PV_Temp, t_NAC_CB_Temp, t_PFC_Pyra_Temp, t_PFC_PV_Temp, t_PFC_CB_Temp};
        Full_Data = {y_Temp};
    end
end
