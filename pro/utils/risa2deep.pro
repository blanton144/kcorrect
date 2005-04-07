; code to convert DEEP filters received from
; Risa (from Chris Willmer originally, I believe
; from Nick Kaiser's web site) to atmospheric 
; version
; assumes 1.3 airmass at Palomar
pro risa2deep

filter=yanny_readone(getenv('KCORRECT_DIR')+ $
                     '/data/filters/risa_deep_B.par', hdr=hdr)
oke, filter.lambda, 1., 1.3, ext
filter.pass=filter.pass*10.^(-0.4*ext)
hdr=[hdr, '# ']
hdr=[hdr, '# Atmospheric extinction added:']
hdr=[hdr, '#    oke, filter.lambda, 1., 1.3, ext']
hdr=[hdr, '#    filter.pass=filter.pass*10.^(-0.4*ext)']
hdr=[hdr, '# MRB, NYU 2005-04-05']
yanny_write, getenv('KCORRECT_DIR')+'/data/filters/deep_B.par', $
  ptr_new(filter), hdr=hdr

filter=yanny_readone(getenv('KCORRECT_DIR')+ $
                     '/data/filters/risa_deep_R.par', hdr=hdr)
oke, filter.lambda, 1., 1.3, ext
filter.pass=filter.pass*10.^(-0.4*ext)
hdr=[hdr, '# ']
hdr=[hdr, '# Atmospheric extinction added:']
hdr=[hdr, '#    oke, filter.lambda, 1., 1.3, ext']
hdr=[hdr, '#    filter.pass=filter.pass*10.^(-0.4*ext)']
hdr=[hdr, '# MRB, NYU 2005-04-05']
yanny_write, getenv('KCORRECT_DIR')+'/data/filters/deep_R.par', $
  ptr_new(filter), hdr=hdr

filter=yanny_readone(getenv('KCORRECT_DIR')+ $
                     '/data/filters/risa_deep_I.par', hdr=hdr)
oke, filter.lambda, 1., 1.3, ext
filter.pass=filter.pass*10.^(-0.4*ext)
hdr=[hdr, '# ']
hdr=[hdr, '# Atmospheric extinction added:']
hdr=[hdr, '#    oke, filter.lambda, 1., 1.3, ext']
hdr=[hdr, '#    filter.pass=filter.pass*10.^(-0.4*ext)']
hdr=[hdr, '# MRB, NYU 2005-04-05']
yanny_write, getenv('KCORRECT_DIR')+'/data/filters/deep_I.par', $
  ptr_new(filter), hdr=hdr

end
