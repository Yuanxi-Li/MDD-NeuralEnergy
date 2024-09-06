function efun = efun_MSN(z)
if (abs(z) < 1e-4)
    efun = 1-z/2;
else
    efun = z/(exp(z)-1);
end
end