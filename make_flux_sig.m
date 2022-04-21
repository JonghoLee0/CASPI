% make an average signal photon vector with a Gaussian shape

function flux_sig = make_flux_sig(photon_sig, bin_N, bin_t, tof_idx, FWHM)

FWHM_bin = FWHM/bin_t;
sigma_bin = FWHM_bin/(2*sqrt(2*log(2)));


% make a Gaussian function with mean and sigma
dummy_x = [1 : bin_N]';
flux_sig = normpdf(dummy_x, tof_idx, sigma_bin);


% make the area equal to signal photons
flux_sig = photon_sig*flux_sig;