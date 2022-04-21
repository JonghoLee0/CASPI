function [SIG_hat, coeff_hard] = shrink_hard_sig(SIG_noisy, th)

coeff_hard = ones(size(SIG_noisy));

idx = (abs(SIG_noisy) <= th);
SIG_noisy(idx) = 0;
coeff_hard(idx) = 0;

SIG_hat = SIG_noisy;