function [k1_BKKCa,k2_BKKCa,k3_BKKCa,k4_BKKCa] = rates_BKKCa_MSN(Vm,Cai)
k1_BKKCa = alp_BKKCa_MSN(0.1,Vm,-10,1);
k2_BKKCa = alp_BKKCa_MSN(0.1,Vm,-120,-10);
k3_BKKCa = alpha_BKKCa_MSN(0.001,1,Vm,-20,7) *1.0e8*(Cai*1.0)^3;
k4_BKKCa = alp_BKKCa_MSN(0.01,Vm,-44,-5);
end