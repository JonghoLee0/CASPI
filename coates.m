function flux_hat = coates(hst, cycle_N)


% if sum(hst) == cycle_N, increment 1 to prevent division by 0
if sum(hst) >= cycle_N
    cycle_N = sum(hst) + 1;
end

cum_sum1 = cumsum(hst);
cum_sum2 = [0; cum_sum1(1:end-1)];

numer = cycle_N - cum_sum2;
denom = cycle_N - cum_sum1;

flux_hat = log(numer./denom);