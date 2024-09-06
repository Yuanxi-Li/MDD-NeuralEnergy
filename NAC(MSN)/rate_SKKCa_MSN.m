function [tau_SKKCa,oinf_SKKCa] = rate_SKKCa_MSN(Vm,Cai);
a_SKKCa = alpha_SKKCa_MSN(Vm,Cai);
b_SKKCa = beta_SKKCa_MSN(Vm,Cai);
tau_SKKCa = 1/(a_SKKCa+b_SKKCa);
oinf_SKKCa = a_SKKCa*tau_SKKCa;
end