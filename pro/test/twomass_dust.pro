pro twomass_dust

savfile='data_twomass_dust.sav'
if(NOT file_test(savfile)) then begin
    im=hogg_mrdfits(vagc_name('object_sdss_imaging'),1,nrowchunk=20000 )
    sp=hogg_mrdfits(vagc_name('object_sdss_spectro'),1,nrowchunk=20000, $
                    columns=['plate', 'fiberid', 'mjd'] )
    kcorrect=hogg_mrdfits(vagc_name('kcorrect',collision_type='none', $
                                    flux_type='petro',band_shift=0.1),1, $
                          nrowchunk=20000)
    twomass=hogg_mrdfits(vagc_name('object_twomass'),1,nrowchunk=20000, $
                         columns=['k_m_k20fe'])
    
; select twomass galaxies
    ii=where(twomass.k_m_k20fe gt 10. and twomass.k_m_k20fe lt 13.5 and $
             kcorrect.z gt 0.01 and kcorrect.z lt 0.1)
    twomass=twomass[ii]
    im=im[ii]
    sp=sp[ii]
    kcorrect=kcorrect[ii]

; read in line info
    readspec, sp.plate, sp.fiberid, mjd=sp.mjd, zline=zline

    calibobj=retrieve_calibobj(im, columns=['ab_exp','ab_dev','fracpsf'])
    save,filename=savfile
endif else begin
    restore,savfile
endelse

ii=where(kcorrect.abmaggies[7] gt 0. and $
         kcorrect.kcorrect[7] eq kcorrect.kcorrect[7])
twomass=twomass[ii]
im=im[ii]
sp=sp[ii]
kcorrect=kcorrect[ii]
calibobj=calibobj[ii]
zline=zline[ii]
dm=lf_distmod(kcorrect.z, omega0=0.3, omegal0=0.7) 
absmk=22.5-2.5*alog10(kcorrect.abmaggies[7])-dm-kcorrect.kcorrect[7]

stop

ii=where(zline[17,*].lineew gt 5.)
absmk=absmk[ii]
twomass=twomass[ii]
im=im[ii]
sp=sp[ii]
kcorrect=kcorrect[ii]


ii=where(absmk gt -22. and absmk lt -21.)
absmk=absmk[ii]
twomass=twomass[ii]
im=im[ii]
sp=sp[ii]
kcorrect=kcorrect[ii]

rmk=-2.5*alog10(kcorrect.abmaggies[2]/kcorrect.abmaggies[7])
umk=-2.5*alog10(kcorrect.abmaggies[0]/kcorrect.abmaggies[7])
zmk=-2.5*alog10(kcorrect.abmaggies[4]/kcorrect.abmaggies[7])

stop

end
