function y = Initial_MSN(V_MSN_Soma,V_MSN_Dendrite);
y = zeros(56,1);
y(1) = V_MSN_Soma;
y(2) = V_MSN_Dendrite;
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


CaL1_2_MSN = MSN_ionmodel2('CaL1.2', 6.7e-6, 0.17, -8.9, -13.4, -6.7, 11.9);
CaL1_3_MSN = MSN_ionmodel2('CaL1.3', 4.25e-7, nan, -33, -13.4, -6.7, 11.9);
CaN_MSN = MSN_ionmodel2('CaN', 1.0e-5, 0.21, -8.7, -74.8, -7.4, 6.5);
CaQ_MSN = MSN_ionmodel2('CaQ', 6.0e-6, nan, -9.0, nan, -6.6, nan);
CaR_MSN = MSN_ionmodel2('CaR', 2.6e-5, nan, -10.3, -33.3, -6.6, 17);
CaT_MSN = MSN_ionmodel2('CaT', 4e-7, nan, -51.73, -80, -6.53, 6.7);

y(3) = infinite(y(1),NaF_MSN_Soma.mVhalf,NaF_MSN_Soma.mk);
y(4) = infinite(y(1),NaF_MSN_Soma.hVhalf,NaF_MSN_Soma.hk);
y(5) = infinite(y(1),NaP_MSN_Soma.mVhalf,NaP_MSN_Soma.mk);
y(6) = infinite(y(1),NaP_MSN_Soma.hVhalf,NaP_MSN_Soma.hk);
y(7) = infinite(y(1),KAf_MSN_Soma.mVhalf,KAf_MSN_Soma.mk);
y(8) = infinite(y(1),KAf_MSN_Soma.hVhalf,KAf_MSN_Soma.hk);
y(9) = infinite(y(1),KAs_MSN_Soma.mVhalf,KAs_MSN_Soma.mk);
y(10) = infinite(y(1),KAs_MSN_Soma.hVhalf,KAs_MSN_Soma.hk);
y(11) = infinite(y(1),KIR_MSN.mVhalf,KIR_MSN.mk);
y(12) = infinite(y(1),KRP_MSN.mVhalf,KRP_MSN.mk);
y(13) = infinite(y(1),KRP_MSN.mVhalf,KRP_MSN.mk);
y(14) = infinite(y(1),CaL1_2_MSN.mVhalf,CaL1_2_MSN.mk);
y(15) = infinite(y(1),CaL1_2_MSN.mVhalf,CaL1_2_MSN.mk);
y(16) = infinite(y(1),CaL1_3_MSN.mVhalf,CaL1_3_MSN.mk);
y(17) = infinite(y(1),CaL1_3_MSN.mVhalf,CaL1_3_MSN.mk);
y(18) = infinite(y(1), CaN_MSN.mVhalf, CaN_MSN.mk);
y(19) = infinite(y(1), CaN_MSN.hVhalf, CaN_MSN.hk);
y(20) = infinite(y(1),CaQ_MSN.mVhalf,CaQ_MSN.mk);
y(21) = infinite(y(1),CaR_MSN.mVhalf,CaR_MSN.mk);
y(22) = infinite(y(1),CaR_MSN.hVhalf,CaR_MSN.hk);
y(23) = infinite(y(1),CaT_MSN.mVhalf,CaT_MSN.mk);
y(24) = infinite(y(1),CaT_MSN.hVhalf,CaT_MSN.hk);
y(25) = 1e-5;%------mmol/L
y(26) = 1e-5;
[tau_SKKCa_MSN_Soma,oinf_SKKCa_MSN_Soma] = rate_SKKCa_MSN(y(1),y(26));
y(27) = oinf_SKKCa_MSN_Soma;
y(28) = 0;
y(29) = 1;
y(30) = 0;
y(31) = infinite(y(2),NaF_MSN_Dendrite.mVhalf,NaF_MSN_Dendrite.mk);
y(32) = infinite(y(2),NaF_MSN_Dendrite.hVhalf,NaF_MSN_Dendrite.hk);
y(33) = infinite(y(2),NaP_MSN_Dendrite.mVhalf,NaP_MSN_Dendrite.mk);
y(34) = infinite(y(2),NaP_MSN_Dendrite.hVhalf,NaP_MSN_Dendrite.hk);
y(35) = infinite(y(2),KAf_MSN_Dendrite.mVhalf,KAf_MSN_Dendrite.mk);
y(36) = infinite(y(2),KAf_MSN_Dendrite.hVhalf,KAf_MSN_Dendrite.hk);
y(37) = infinite(y(2),KAs_MSN_Dendrite.mVhalf,KAs_MSN_Dendrite.mk);
y(38) = infinite(y(2),KAs_MSN_Dendrite.hVhalf,KAs_MSN_Dendrite.hk);
y(39) = infinite(y(2),KIR_MSN.mVhalf,KIR_MSN.mk);
y(40) = infinite(y(2),CaL1_2_MSN.mVhalf,CaL1_2_MSN.mk);
y(41) = infinite(y(2),CaL1_2_MSN.mVhalf,CaL1_2_MSN.mk);
y(42) = infinite(y(2),CaL1_3_MSN.mVhalf,CaL1_3_MSN.mk);
y(43) = infinite(y(2),CaL1_3_MSN.mVhalf,CaL1_3_MSN.mk);
y(44) = infinite(y(2), CaN_MSN.mVhalf, CaN_MSN.mk);
y(45) = infinite(y(2), CaN_MSN.hVhalf, CaN_MSN.hk);
y(46) = infinite(y(2),CaQ_MSN.mVhalf,CaQ_MSN.mk);
y(47) = infinite(y(2),CaR_MSN.mVhalf,CaR_MSN.mk);
y(48) = infinite(y(2),CaR_MSN.hVhalf,CaR_MSN.hk);
y(49) = infinite(y(2),CaT_MSN.mVhalf,CaT_MSN.mk);
y(50) = infinite(y(2),CaT_MSN.hVhalf,CaT_MSN.hk);
y(51) = 1e-5;%------mmol/L
y(52) = 1e-5;
[tau_SKKCa_MSN_Dendrite,oinf_SKKCa_MSN_Dendrite] = rate_SKKCa_MSN(y(2),y(52));
y(53) = oinf_SKKCa_MSN_Dendrite;
y(54) = 0;
y(55) = 1;
y(56) = 0;









