pro andrey_dust

savfile='data_andrey_dust.sav'
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
    ;readspec, sp.plate, sp.fiberid, mjd=sp.mjd, zline=zline

    ;calibobj=retrieve_calibobj(im, columns=['ab_exp','ab_dev','fracpsf'])
    save,filename=savfile
endif else begin
    restore,savfile
endelse

ii=where((im.vagc_select and 4) gt 0 and $
         kcorrect.abmaggies[7] gt 0. and $
         kcorrect.kcorrect[7] eq kcorrect.kcorrect[7])
twomass=twomass[ii]
im=im[ii]
sp=sp[ii]
kcorrect=kcorrect[ii]
dm=lf_distmod(kcorrect.z, omega0=0.3, omegal0=0.7) 
absmk=22.5-2.5*alog10(kcorrect.abmaggies[7])-dm-kcorrect.kcorrect[7]
rmk=-2.5*alog10(kcorrect.abmaggies[2]/kcorrect.abmaggies[7])- $
  kcorrect.kcorrect[2]+kcorrect.kcorrect[7]
gmr=-2.5*alog10(kcorrect.abmaggies[1]/kcorrect.abmaggies[2])- $
  kcorrect.kcorrect[1]+kcorrect.kcorrect[2]

ndata=n_elements(absmk)
ndims=3
data=fltarr(ndims,ndata)
data[0,*]=absmk
data[1,*]=rmk
data[2,*]=gmr
range=[[-17.,-23.0], $
       [-0.8,2.5], $
       [0.1,1.2]]
label=['!8M!dK!n!6', $
       '!8!u0.1!n[r-K]!6', $
       '!8!u0.1!n[g-r]!6' ]

jj=lindgen(ndata)
hogg_manyd_scatterplot, fltarr(ndata)+1., data[*,jj], 'andrey_dust.ps', $
  range=range,/conditional,label=label, xnpix=30, ynpix=30, $
  quantiles=[0.025, 0.16, 0.5, 0.84, 0.975], satfrac=0.03, manyd=manyd
mwrfits,reform(manyd[*,*,0,1],30,30),'andrey_dust.fits',/create
k_write_ascii_table,reform(manyd[*,*,0,1],30,30),'andrey_dust.dat'
openw,unit,'andrey_dust_data.dat',/get_lun
for i=0L, ndata-1L do $
  printf,unit,strtrim(string(absmk[i],format='(f40.5)'),2)+' '+ $
  strtrim(string(rmk[i],format='(f40.5)'),2)
free_lun,unit

stop

end
