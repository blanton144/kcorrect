pro twomass_dust

im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1,nrowchunk=20000 )
sp=hogg_mrdfits(vagc_name('object_sdss_spectro'),1,nrowchunk=20000, $
               columns=['plate', 'fiberid', 'mjd'] )
kcorrect=hogg_mrdfits(vagc_name('kcorrect',collision_type='none', $
                                flux_type='petro',band_shift=0.1),1, $
                      nrowchunk=20000)
twomass=hogg_mrdfits(vagc_name('object_twomass'),1,nrowchunk=20000, $
                     columns=['k_m_k20fe'])

ii=where(twomass.k_m_k20fe gt 10. and twomass.k_m_k20fe lt 13.5 and $
         kcorrect.z gt 0.01 and kcorrect.z lt 0.1)
twomass=twomass[ii]
im=im[ii]
sp=sp[ii]
kcorrect=kcorrect[ii]

dm=lf_distmod(kcorrect.z, omega0=0.3, omegal0=0.7) 
absmk=22.5-2.5*alog10(kcorrect.abmaggies[7])-dm-kcorrect.kcorrect[7]

ii=where(absmk gt -22. and absmk lt -21.)
absmk=absmk[ii]
twomass=twomass[ii]
im=im[ii]
sp=sp[ii]
kcorrect=kcorrect[ii]

readspec, sp.plate, sp.fiberid, mjd=sp.mjd, zline=zline
ii=where(zline[17,*].lineew gt 5.)
absmk=absmk[ii]
twomass=twomass[ii]
im=im[ii]
sp=sp[ii]
kcorrect=kcorrect[ii]

calibobj=retrieve_calibobj(im, columns='ab_exp')
rmk=-2.5*alog10(kcorrect.abmaggies[2]/kcorrect.abmaggies[7])
umk=-2.5*alog10(kcorrect.abmaggies[0]/kcorrect.abmaggies[7])
zmk=-2.5*alog10(kcorrect.abmaggies[4]/kcorrect.abmaggies[7])

stop

end
