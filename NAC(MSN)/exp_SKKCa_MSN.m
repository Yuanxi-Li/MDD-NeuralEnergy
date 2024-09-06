function exp_SKKCa = exp_SKKCa_MSN(k,d,Vm);
F = 96489;%-------------Faraday(C/mod)
R = 8.31;%------------J/(mol*K)
T = 35+273.15;
exp_SKKCa = k*exp(-2*d*F*Vm/R/T);
end