function flux_hat_grp = coates_grp(hst_grp, N_cycle)

x_size = size(hst_grp, 2);
y_size = size(hst_grp, 1);

flux_hat_grp = zeros(size(hst_grp));

for y = 1 : y_size
    for x = 1 : x_size

        hst = squeeze(hst_grp(y, x, :));
        
        flux_hat = coates(hst, N_cycle);
        flux_hat(flux_hat<0) = 0;
        flux_hat_grp(y, x, :) = flux_hat;

    end
end