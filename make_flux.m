function flux = make_flux(photon_sig, photon_amb, bin_N, bin_t, tof_idx, FWHM)


% average signal photon number vector
flux_sig = make_flux_sig(photon_sig, bin_N, bin_t, tof_idx, FWHM);


% average ambient photon number vector
flux_amb = photon_amb*ones(bin_N, 1);


% average total photon number vector
flux = flux_sig + flux_amb;