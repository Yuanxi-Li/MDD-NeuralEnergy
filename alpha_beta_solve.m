function [mhinf, mhtau] = alpha_beta_solve(mhalpha, mhbeta)
% use m/h alpha and beta to calculate inf and tau
mhinf = mhalpha / (mhalpha + mhbeta);
mhtau = 1 / (mhalpha + mhbeta);
end