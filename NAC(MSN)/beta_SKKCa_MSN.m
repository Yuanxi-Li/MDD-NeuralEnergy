function beta_SKKCa = beta_SKKCa_MSN(Vm,Cai);
k2_SKKCa = 0.011;
d2_SKKCa = 1.0;
bbar_SKKCa = 0.28;
beta_SKKCa = bbar_SKKCa/(1 + Cai/exp_SKKCa_MSN(k2_SKKCa,d2_SKKCa,Vm));
end