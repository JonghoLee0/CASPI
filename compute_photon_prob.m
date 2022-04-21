function photon_prob = compute_photon_prob(avg_photon_vect)

sumcum = cumsum(avg_photon_vect);
sumcum = [0; sumcum(1:end-1)];

photon_prob = (1 - exp(-avg_photon_vect)).*exp(-sumcum);