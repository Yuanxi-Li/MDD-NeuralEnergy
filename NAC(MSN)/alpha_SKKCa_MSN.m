function alpha_SKKCa = alpha_SKKCa_MSN(Vm,Cai);
k1_SKKCa = 0.18;
d1_SKKCa = 0.84;
abar_SKKCa = 0.48;
alpha_SKKCa = abar_SKKCa/(1 + exp_SKKCa_MSN(k1_SKKCa,d1_SKKCa,Vm)/Cai);
end