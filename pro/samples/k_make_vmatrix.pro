pro k_make_vmatrix

k_make_sdss_training_set
k_make_templates

for i=2, 5 do $
  k_fit_sdss_training_set,navg=i

end
