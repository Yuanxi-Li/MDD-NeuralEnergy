function g_Ca = g_Ca_MSN(Vm,ci,co);
F = 96489;%-------------Faraday(C/mod)
R = 8.31;%------------J/(mol*K)
T = 35+273.15;%--------K
z = (1e-3)*2*F*Vm/(R*T);
eco = co*efun_MSN(z);
eci = ci*efun_MSN(-z);
g_Ca = (0.001)*2*F*(eci-eco);
end