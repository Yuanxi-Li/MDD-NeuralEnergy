function current_MSN_record
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

I_NaF_MSN_Soma = NaF_MSN_Soma.gmax.*(y(:,3).^3).*y(:,4).*(y(:,1)-NaF_MSN_Soma.rev);%-----mA
I_NaP_MSN_Soma = NaP_MSN_Soma.gmax.*y(:,5).*y(:,6).*(y(:,1)-NaP_MSN_Soma.rev);%---------mA
I_KAf_MSN_Soma = KAf_MSN_Soma.gmax.*(y(:,7).^2).*y(:,8).*(y(:,1)-KAf_MSN_Soma.rev);%---------mA
I_KAs_MSN_Soma = KAs_MSN_Soma.gmax.*(y(:,9).^2).*(KAs_MSN_Soma.a.*y(10)+(1-KAs_MSN_Soma.a)).*(y(1)-KAs_MSN_Soma.rev);%---------mA
I_KIR_MSN_Soma = KIR_MSN.gmax.*y(:,11).*(y(:,1)-KIR_MSN.rev);%---------mA
I_KRP_MSN_Soma = KRP_MSN.gmax.*y(:,12).*(KRP_MSN.a.*y(:,13)+(1-KRP_MSN.a)).*(y(:,1)-KRP_MSN.rev);%---------mA
Ca_0_MSN = 5;
CaL_0_MSN = 5;
for k = 1:length(y(:,1))
g_CaL1_2_MSN_Soma(k) = g_Ca_MSN(y(k,1),y(k,25),CaL_0_MSN);
g_CaL1_3_MSN_Soma(k) = g_Ca_MSN(y(k,1),y(k,25),CaL_0_MSN);
g_CaN_MSN_Soma(k) = g_Ca_MSN(y(k,1),y(k,26),Ca_0_MSN);
g_CaQ_MSN_Soma(k) = g_Ca_MSN(y(k,1),y(k,26),Ca_0_MSN);
g_CaR_MSN_Soma(k) = g_Ca_MSN(y(k,1),y(k,26),Ca_0_MSN);
g_CaT_MSN_Soma(k) = g_Ca_MSN(y(k,1),y(k,25),Ca_0_MSN);

end
I_CaL1_2_MSN_Soma = g_CaL1_2_MSN_Soma'.*CaL1_2_MSN.Pbar.*y(:,14).*y(:,14).*(y(:,15).*CaL1_2_MSN.a+(1-CaL1_2_MSN.a));%------mA

I_CaL1_3_MSN_Soma = g_CaL1_3_MSN_Soma'.*CaL1_3_MSN.Pbar.*y(:,16).*y(:,16).*y(:,17);%------mA
I_CaN_MSN_Soma = g_CaN_MSN_Soma'.*CaN_MSN.Pbar.*y(:,18).*y(:,18).*(CaN_MSN.a.*y(:,19)+(1-CaN_MSN.a));%------mA
I_CaQ_MSN_Soma = g_CaQ_MSN_Soma'.*CaQ_MSN.Pbar.*y(:,20).*y(:,20);%------mA
I_CaR_MSN_Soma = g_CaR_MSN_Soma'.*CaR_MSN.Pbar.*y(:,21).*y(:,21).*y(:,21).*y(:,22);%------mA
I_CaT_MSN_Soma = g_CaT_MSN_Soma'.*CaT_MSN.Pbar.*y(:,23).*y(:,23).*y(:,23).*y(:,24);%------mA
I_SKKCa_MSN_Soma = SKKCa_MSN.gmax.*y(:,27).*(y(:,1)-SKKCa_MSN.rev);
I_BKKCa_MSN_Soma = BKKCa_MSN.gmax.*y(:,28).*(y(:,1)-BKKCa_MSN.rev);





