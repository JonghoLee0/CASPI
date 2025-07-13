function SIG_cubelet_hat = shrink_hard(SIG_cubelet, th)

idx = (abs(SIG_cubelet) <= th);
SIG_cubelet(idx) = 0;

SIG_cubelet_hat = SIG_cubelet;