function [t,y,t_NAC_MSN,t_NAC_PV,t_NAC_CB,t_PFC_Pyra,t_PFC_PV,t_PFC_CB] = FiringSequence_Dopamineratio_Change(Seed_Net, Seed_Ext, Dopamine_Ratio, I_Ext_MSN_Soma, I_Ext_MSN_Dendrite,...
    I_Ext_PFC_Pyra_Soma,I_Ext_PFC_Pyra_Proximal,I_Ext_PFC_Pyra_Distal,...
    I_Ext_PFC_Inter_CB,I_Ext_PFC_Inter_PV, I_Ext_NAc_Inter_CB, I_Ext_NAc_Inter_PV)

[t,y] = NAC_PFC_Dopamine_Ratio_Change(Seed_Net, Seed_Ext, Dopamine_Ratio, I_Ext_MSN_Soma, I_Ext_MSN_Dendrite,...
    I_Ext_PFC_Pyra_Soma,I_Ext_PFC_Pyra_Proximal,I_Ext_PFC_Pyra_Distal,...
    I_Ext_PFC_Inter_CB,I_Ext_PFC_Inter_PV, I_Ext_NAc_Inter_CB, I_Ext_NAc_Inter_PV);


y1 = y;
[m_NAC_MSN,n_NAC_MSN] = findpeaks(y1(:,1));
n_NAC_MSN = n_NAC_MSN(m_NAC_MSN>0);
peaks_NAC_MSN = y1(n_NAC_MSN,1);
t_NAC_MSN = t(n_NAC_MSN);
[m_NAC_PV,n_NAC_PV] = findpeaks(y1(:,99));
n_NAC_PV = n_NAC_PV(m_NAC_PV>0);
peaks_NAC_PV = y1(n_NAC_PV,99);
t_NAC_PV = t(n_NAC_PV);
[m_NAC_CB,n_NAC_CB] = findpeaks(y1(:,143));
n_NAC_CB = n_NAC_CB(m_NAC_CB>0);
peaks_NAC_CB = y1(n_NAC_CB,143);
t_NAC_CB = t(n_NAC_CB);
for k = 0:19
[m_PFC_Pyra{k+1},n_PFC_Pyra{k+1}] = findpeaks(y1(:,189+55*k));
n_PFC_Pyra{k+1} = n_PFC_Pyra{k+1}(m_PFC_Pyra{k+1}>0);
peaks_PFC_Pyra{k+1} = y1(n_PFC_Pyra{k+1},189+55*k);
t_PFC_Pyra{k+1} = t(n_PFC_Pyra{k+1});
end
for k = 0:2
[m_PFC_PV{k+1},n_PFC_PV{k+1}] = findpeaks(y1(:,1289+47*k));
n_PFC_PV{k+1} = n_PFC_PV{k+1}(m_PFC_PV{k+1}>0);
peaks_PFC_PV{k+1} = y1(n_PFC_PV{k+1},1289+47*k);
t_PFC_PV{k+1} = t(n_PFC_PV{k+1});
end
for k = 0:1
[m_PFC_CB{k+1},n_PFC_CB{k+1}] = findpeaks(y1(:,1430+49*k));
n_PFC_CB{k+1} = n_PFC_CB{k+1}(m_PFC_CB{k+1}>0);
peaks_PFC_CB{k+1} = y1(n_PFC_CB{k+1},1430+49*k);
t_PFC_CB{k+1} = t(n_PFC_CB{k+1});
end
% figure()
% size = 10;
% scatter(t_NAC_MSN,ones(length(t_NAC_MSN),1),'filled','MarkerFaceColor',[255 103 102]/255,'SizeData',size)
% hold on
% scatter(t_NAC_PV,2*ones(length(t_NAC_PV),1),'filled','MarkerFaceColor',[52 102 103]/255,'SizeData',size)
% scatter(t_NAC_CB,3*ones(length(t_NAC_CB),1),'filled','MarkerFaceColor',[155 153 104]/255,'SizeData',size)
% for k = 0:19
% scatter(t_PFC_Pyra{k+1},(k+4)*ones(length(t_PFC_Pyra{k+1}),1),'filled','MarkerFaceColor',[131 101 106]/255,'SizeData',size)
% end
% for k = 0:2
% scatter(t_PFC_PV{k+1},(k+24)*ones(length(t_PFC_PV{k+1}),1),'filled','MarkerFaceColor',[102 155 205]/255,'SizeData',size)
% end
% for k = 0:1
% scatter(t_PFC_CB{k+1},(k+27)*ones(length(t_PFC_CB{k+1}),1),'filled','MarkerFaceColor',[103 103 155]/255,'SizeData',size)
% end

% size = 10;
% scatter(Full_Range{1,1}{1,1},ones(length(Full_Range{1,1}{1,1}),1),'filled','MarkerFaceColor',[255 103 102]/255,'SizeData',size)
% hold on
% scatter(Full_Range{1,2}{1,1},2*ones(length(Full_Range{1,2}{1,1}),1),'filled','MarkerFaceColor',[52 102 103]/255,'SizeData',size)
% scatter(Full_Range{1,3}{1,1},3*ones(length(Full_Range{1,3}{1,1}),1),'filled','MarkerFaceColor',[155 153 104]/255,'SizeData',size)
% for k = 0:19
% scatter(Full_Range{1,4}{1,1}{k+1},(k+4)*ones(length(Full_Range{1,4}{1,1}{k+1}),1),'filled','MarkerFaceColor',[131 101 106]/255,'SizeData',size)
% end
% for k = 0:2
% scatter(Full_Range{1,5}{1,1}{k+1},(k+24)*ones(length(Full_Range{1,5}{1,1}{k+1}),1),'filled','MarkerFaceColor',[102 155 205]/255,'SizeData',size)
% end
% for k = 0:1
% scatter(Full_Range{1,6}{1,1}{k+1},(k+27)*ones(length(Full_Range{1,6}{1,1}{k+1}),1),'filled','MarkerFaceColor',[103 103 155]/255,'SizeData',size)
% end