function alp_BKKCa = alp_BKKCa_MSN(tmin,Vm,Vhalf,k)
alp_BKKCa = 1/(tmin+exp((Vm-Vhalf)/k)*1.0);
end