function mhinf = infinite(Vm, Vhalf, k);
% minf = 1/(1+exp((V-mVhalf)/mk)),
% hinf = 1/(1+exp((V-hVhalf)/hk))
mhinf = 1/(1+exp((Vm-Vhalf)/k));
end