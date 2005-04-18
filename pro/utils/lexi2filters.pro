pro lexi2filters

restore, getenv('KCORRECT_DIR')+'/data/filters/goods/goodsfilters.sav'

hdr0='# Units:'
hdr0=[hdr0, '#  "lambda" is in Angstroms']
hdr0=[hdr0, '#  "pass" is the contribution to the detector signal per photon']
hdr0=[hdr0, '#']
hdr0=[hdr0, '# Bandpass Name(s): ']
hdr0=[hdr0, '#']
hdr0=[hdr0, '# Instrument: ISAAC']
hdr0=[hdr0, '#']
hdr0=[hdr0, '# Determined by: Lexi Moustakas?']
hdr0=[hdr0, '#']
hdr0=[hdr0, '# Date of determination: ???']
hdr0=[hdr0, '#']
hdr0=[hdr0, '# Notes:']
hdr0=[hdr0, '#']
hdr0=[hdr0, '# From Lexi Moustakas -- he will get me provenance soon?']
hdr0=[hdr0, '#  (2005-04-18)']

kfilter1={lambda:0.D, $
          pass:0.D}

hdr=hdr0
ii=where(hdr eq '# Bandpass Name(s): ')
hdr[ii]='# Bandpass Names(s): H'
kfilter=replicate(kfilter1, n_elements(ih.filtw))
kfilter.lambda=ih.filtw
kfilter.pass=ih.filtf
yanny_write, getenv('KCORRECT_DIR')+'/data/filters/H_ISAAC.par', $
  ptr_new(kfilter), hdr=hdr

hdr=hdr0
ii=where(hdr eq '# Bandpass Name(s): ')
hdr[ii]='# Bandpass Names(s): J'
kfilter=replicate(kfilter1, n_elements(ij.filtw))
kfilter.lambda=ij.filtw
kfilter.pass=ij.filtf
yanny_write, getenv('KCORRECT_DIR')+'/data/filters/J_ISAAC.par', $
  ptr_new(kfilter), hdr=hdr

hdr=hdr0
ii=where(hdr eq '# Bandpass Name(s): ')
hdr[ii]='# Bandpass Names(s): Ks'
kfilter=replicate(kfilter1, n_elements(iks.filtw))
kfilter.lambda=iks.filtw
kfilter.pass=iks.filtf
yanny_write, getenv('KCORRECT_DIR')+'/data/filters/Ks_ISAAC.par', $
  ptr_new(kfilter), hdr=hdr

hdr0='# Units:'
hdr0=[hdr0, '#  "lambda" is in Angstroms']
hdr0=[hdr0, '#  "pass" is the contribution to the detector signal per photon']
hdr0=[hdr0, '#']
hdr0=[hdr0, '# Bandpass Name(s): ']
hdr0=[hdr0, '#']
hdr0=[hdr0, '# Instrument: HST ACS']
hdr0=[hdr0, '#']
hdr0=[hdr0, '# Determined by: ???']
hdr0=[hdr0, '#']
hdr0=[hdr0, '# Date of determination: ???']
hdr0=[hdr0, '#']
hdr0=[hdr0, '# Notes:']
hdr0=[hdr0, '#']
hdr0=[hdr0, '# From Lexi Moustakas -- he will get me provenance soon?']
hdr0=[hdr0, '#  (2005-04-18)']

hdr=hdr0
ii=where(hdr eq '# Bandpass Name(s): ')
hdr[ii]='# Bandpass Names(s): F435W (B)'
kfilter=replicate(kfilter1, n_elements(ab.filtw))
kfilter.lambda=ab.filtw
kfilter.pass=ab.filtf
yanny_write, getenv('KCORRECT_DIR')+'/data/filters/acs_f435w.par', $
  ptr_new(kfilter), hdr=hdr

hdr=hdr0
ii=where(hdr eq '# Bandpass Name(s): ')
hdr[ii]='# Bandpass Names(s): F775W (i)'
kfilter=replicate(kfilter1, n_elements(ai.filtw))
kfilter.lambda=ai.filtw
kfilter.pass=ai.filtf
yanny_write, getenv('KCORRECT_DIR')+'/data/filters/acs_f775w.par', $
  ptr_new(kfilter), hdr=hdr

hdr=hdr0
ii=where(hdr eq '# Bandpass Name(s): ')
hdr[ii]='# Bandpass Names(s): F606W (V)'
kfilter=replicate(kfilter1, n_elements(av.filtw))
kfilter.lambda=av.filtw
kfilter.pass=av.filtf
yanny_write, getenv('KCORRECT_DIR')+'/data/filters/acs_f606w.par', $
  ptr_new(kfilter), hdr=hdr

hdr=hdr0
ii=where(hdr eq '# Bandpass Name(s): ')
hdr[ii]='# Bandpass Names(s): F850lp (z)'
kfilter=replicate(kfilter1, n_elements(az.filtw))
kfilter.lambda=az.filtw
kfilter.pass=az.filtf
yanny_write, getenv('KCORRECT_DIR')+'/data/filters/acs_f850lp.par', $
  ptr_new(kfilter), hdr=hdr

end
