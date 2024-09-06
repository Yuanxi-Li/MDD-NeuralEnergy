%% Initialization (Adjust the parameters to get different results)

% Set up stimuli currents 
I_Ext_MSN_Soma= ? ;
I_Ext_MSN_Dendrite= ? ; 
I_Ext_PFC_Pyra_Soma= ? ; % Here, I_Pyra_Proximal and I_Pyra_Distal was set equally to I_Pyra_Soma
I_Ext_PFC_Inter_CB= ? ; 
I_Ext_PFC_Inter_PV= ? ; 
I_Ext_NAc_Inter_CB= ? ; 
I_Ext_NAc_Inter_PV= ? ; 

% Randomization Parameters
Seed_Ext = ? ;
Seed_Net = ? ;

% Set up dopamine concentration
Seed_Dopamine = ? ;
rng(Seed_Dopamine);
Dopamine_Ratio_Range = ? ; 
Dopamine_Ratio = randi([Dopamine_Ratio_Range*1e6 (Dopamine_Ratio_Range+0.25)*1e6],1,3005/0.02)/1e6; % It's used to generate a random time-series vector to mimic dopamine concentration.

% Simulation
[t,y,t_NAC_MSN,t_NAC_PV,t_NAC_CB,t_PFC_Pyra,t_PFC_PV,t_PFC_CB] = FiringSequence_Dopamineratio_Change(Seed_Net, Seed_Ext, Dopamine_Ratio,...
I_Ext_MSN_Soma,I_Ext_MSN_Dendrite,I_Ext_PFC_Pyra_Soma,I_Ext_PFC_Pyra_Soma,I_Ext_PFC_Pyra_Soma,...
I_Ext_PFC_Inter_CB,I_Ext_PFC_Inter_PV, I_Ext_NAc_Inter_CB, I_Ext_NAc_Inter_PV);