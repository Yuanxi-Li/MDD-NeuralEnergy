function y0 = Initial_NAc_PFC(V_PFC_Inter_PV,V_PFC_Pyra_Soma,V_PFC_Pyra_Proximal,V_PFC_Pyra_Distal,...
    V_MSN_Soma,V_MSN_Dendrite,V_PFC_Inter_CB,Dopamine_Ratio)
y0 = zeros(1527,1);
y1 = Initial_MSN(V_MSN_Soma,V_MSN_Dendrite);
y2 = Initial_PFC_Inter_PV(V_PFC_Inter_PV);
y3 = Initial_PFC_Inter_CB(V_PFC_Inter_CB);
y4 = Initial_PFC_Pyra(V_PFC_Pyra_Soma,V_PFC_Pyra_Proximal,V_PFC_Pyra_Distal,Dopamine_Ratio);

    

y0(1:56) = y1;
y0(99:101) = y2;
y0(143:147) = y3;

for Num_PFC_Pyra = 0:19
    y0(189+55*Num_PFC_Pyra:200+55*Num_PFC_Pyra) = y4;
end
for Num_PFC_PV = 0:2
    y0(1289+Num_PFC_PV*47:1291+Num_PFC_PV*47) = y2;
end
for Num_PFC_CB = 0:1
    y0(1430+49*Num_PFC_CB:1434+49*Num_PFC_CB) = y3;
end






end