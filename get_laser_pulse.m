function psf = get_laser_pulse(sigma_t, delta_t)

sigma_bin = sigma_t/delta_t;                
bin_max = 10*round(sigma_bin);
psf = normpdf([1:bin_max]', bin_max/2, sigma_bin);