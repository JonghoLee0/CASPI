function flux_hat = coates(hst, N_cycle)


% if sum(hst) == cycle_N, increment 1 to prevent division by 0
if sum(hst) >= N_cycle
    N_cycle = sum(hst) + 1;
end

cum_sum1 = cumsum(hst);
cum_sum2 = [0; cum_sum1(1:end-1)];

numer = N_cycle - cum_sum2;
denom = N_cycle - cum_sum1;

flux_hat = log(numer./denom);