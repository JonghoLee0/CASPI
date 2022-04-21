function [hst, norm_hst, photon_prob] = make_hst(flux, cycle_N)


% identify number of bins
bin_N = size(flux, 1);


% create noisy photon vector set
noisy_flux_set = poissrnd(repmat(flux, 1, cycle_N));


% find the first photon index
[max_val, max_idx] = max(noisy_flux_set~=0, [], 1);
first_photon_idx = max_idx(max_val==1);



% create histogram with the first photon indices
hst = histcounts(first_photon_idx, [0.5 : 1 : bin_N+0.5]);
hst = hst';



% normalized histogram
norm_hst = hst/cycle_N;



% compute photon detection probability
photon_prob = compute_photon_prob(flux);