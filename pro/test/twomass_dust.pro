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
zline=zline[*,ii]
dm=lf_distmod(kcorrect.z, omega0=0.3, omegal0=0.7) 
absmk=22.5-2.5*alog10(kcorrect.abmaggies[7])-dm-kcorrect.kcorrect[7]
rmk=-2.5*alog10(kcorrect.abmaggies[2]/kcorrect.abmaggies[7])- $
  kcorrect.kcorrect[2]+kcorrect.kcorrect[7]
gmr=-2.5*alog10(kcorrect.abmaggies[1]/kcorrect.abmaggies[2])- $
  kcorrect.kcorrect[1]+kcorrect.kcorrect[2]

ndata=n_elements(absmk)
ndims=5
data=fltarr(ndims,ndata)
data[0,*]=calibobj.ab_exp[2]
data[1,*]=calibobj.fracpsf[2]
data[2,*]=absmk
data[3,*]=zline[17,*].lineew
data[4,*]=rmk
range=[[0.,1.], $
       [0.,1.], $
       [-23.,-17.], $
       [-2.,13.], $
       [0.5,2.]]
labels=['!8b/a!6', $
        '!8f_{!6dev}', $
        '!8M_{K}!6', $
        '!6H_{\alpha}!6', $
        '!6H_{\alpha}/H_{\beta}!6', $
        '!8!u0.1!n[r-K]!6']

stop

jj=lindgen(ndata)
hogg_manyd_scatterplot, fltarr(ndata)+1., data[*,jj], 'color.ps', range=range
hogg_manyd_meanplot, fltarr(ndata)+1., data[*,jj], ndims-1L, $
  'color_vs_ba.ps', $
  xdims=lindgen(ndims-1),ydims=lindgen(ndims-1), levels=0.1*findgen(20), $
  dbin=[0.2,0.1,1.0,2.,0.1], range=range

nn=6
for i=1L, 3L do begin
    for j=0L, nn-1L do begin
        lo=range[0,i]+(range[1,i]-range[0,i])*float(j)/float(nn)
        hi=range[0,i]+(range[1,i]-range[0,i])*float(j+1)/float(nn)
        jj=where(data[i,*] gt lo and data[i,*] lt hi, jjcount)
        if(jjcount gt 0) then begin
            ndata=n_elements(jj)
            hogg_manyd_meanplot, fltarr(ndata)+1., data[*,jj], ndims-1L, $
              'color_vs_ba_'+strtrim(string(i),2)+'_'+ $
              strtrim(string(j),2)+'.ps', $
              xdims=lindgen(ndims-1),ydims=lindgen(ndims-1), $
              levels=0.1*findgen(20), $
              dbin=[0.2,0.1,1.0,3.,0.1], range=range
        endif
    endfor
endfor

stop

ii=lindgen(n_elements(absmk))
plot,calibobj[ii].ab_exp[2],calibobj[ii].ab_dev[2],psym=4

ii=where(zline[17,*].lineew gt 5.)
ii=where(calibobj.fracpsf[2] lt 0.5)
rmk=-2.5*alog10(kcorrect[ii].abmaggies[2]/kcorrect[ii].abmaggies[7])
plot,calibobj[ii].ab_exp[2],rmk,psym=3

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

red_fac = [5.155, 3.793, 2.751, 2.086, 1.479 ]
red_fac=[red_fac,  0.902, 0.576, 0.367]
ext_cont=fltarr(8,n_elements(spec))
for i=0, 7 do ext_cont[i,*]=spec.tauv_cont*red_fac[i]*2.5/alog(10.)
corrmaggies=kcorrect[m2].abmaggies
for i=0, 7 do corrmaggies[i,*]=corrmaggies[i,*]*10.^(0.4*ext_cont[i,m1])

corrrmk=-2.5*alog10(corrmaggies[2]/corrmaggies[7])

stop

end
