function alpha_BKKCa = alpha_BKKCa_MSN(tmin,tmax,Vm,Vhalf,k)
alpha_BKKCa = 1/(tmin+1/(1/(tmax-tmin)+exp((Vm-Vhalf)/k)*1));
end