function [g_Rall_MSN] = MSN_RallModel
Soma_MSN = MSNComponent('Soma', 16, 16, 1);
Proximal_MSN = MSNComponent('Proximal', 20, 2.25, 1);

% r_L: axial resistance, ¦¸/cm
r_L_MSN = 100;
%(1):SomaByProximal,(2):ProximalBySoma,
%(3):ProximalByMiddle,(4):MiddleByProximal,
%(5):MiddleByDistal,(6):DistalByMiddle
g_Rall_MSN(1) = 10^4 * (Soma_MSN.Diameter * Proximal_MSN.Diameter^2) / ...
(r_L_MSN * Soma_MSN.Length * (Soma_MSN.Diameter^2 * Proximal_MSN.Length + Proximal_MSN.Diameter^2 * Soma_MSN.Length));
g_Rall_MSN(2) = 10^4 * (Proximal_MSN.Diameter * Soma_MSN.Diameter^2) / ...
(r_L_MSN * Proximal_MSN.Length * (Proximal_MSN.Diameter^2 * Soma_MSN.Length + Soma_MSN.Diameter^2 * Proximal_MSN.Length));
end