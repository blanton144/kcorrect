; make demos for web page
pro kc2vl, spec, vmatrix, lambda
vmatrix=spec.flux
lambda=fltarr(n_elements(vmatrix)+1L)
dlambda=spec[1].wavelength-spec[0].wavelength
lambda[0:n_elements(vmatrix)-1]=spec.wavelength-0.5*dlambda
lambda[n_elements(vmatrix)]=spec[n_elements(vmatrix)-1L].wavelength+0.5*dlambda
end
;
pro mk_demo_kcorrect

templates=['elliptical_template.fits', $
           's0_template.fits', $
           'sa_template.fits', $
           'sb_template.fits', $
           'starb3_template.fits'] 

filterlist=['sdss_u0.par','sdss_g0.par','sdss_r0.par','sdss_i0.par']
for i=0L, n_elements(templates)-1L do begin
    spec=mrdfits(getenv('KCORRECT_DIR')+'/data/seds/kc96/'+templates[i],1)
    kc2vl,spec,vmatrix,lambda
    zz=0.25*(findgen(100)+0.5)/100.
    coeffs=replicate(1.,100)
    k_reconstruct_maggies, coeffs, zz, recmaggies, $
      vmatrix=vmatrix, lambda=lambda,filterlist=filterlist

    kcorrect,dum1,dum2,zz,coeffs=coeffs,vmatrix=vmatrix,lambda=lambda, $
      tmpkcorrect1,filterlist=filterlist, band_shift=0.1
    
    kcorrect,recmaggies,recmaggies*0.1,zz,tmpkcorrect2,coeffs=coeffs, $
      /stddev,filterlist=filterlist,vmatrix=vv, lambda=ll, band_shift=0.1, $
      vfile='vmatrix.smooth3.dat', lfile='lambda.smooth3.dat', chi2=chi2
    k_reconstruct_maggies, coeffs, zz, recmaggies2, $
      vmatrix=vv, lambda=ll,filterlist=filterlist
    stop
    
endfor

stop

plates=455+lindgen(1)
fiberids=lindgen(640)+1L
gals1={plate:0L, fiberid:0L}
for i=0L, n_elements(plates)-1L do begin
    readspec,plates[i],fiberids,zans=zans
    ii=where(zans.z gt 0.005 and zans.z lt 0.7 and zans.zwarning eq 0 and $
             zans.class eq 'GALAXY', iicount) 
    if(iicount gt 0) then begin
        tmpgals=replicate(gals1,iicount)
        tmpgals.plate=plates[i]
        tmpgals.fiberid=fiberids[ii]
        if(n_tags(gals) eq 0) then $
          gals=tmpgals $
        else $
          gals=[gals,tmpgals]
    endif
endfor

kcorrect=fltarr(3,n_elements(gals))
kcorrect_real=fltarr(3,n_elements(gals))
for ifiberid=0, n_elements(gals)-1L do begin
    plate=gals[ifiberid].plate
    fiberid=gals[ifiberid].fiberid
    readspec,plate,fiberid,zans=zans,flux=flux,wave=wave
    flux=smooth(flux,20)
    kcorrect,zans.counts_spectro[1:3],zans.counts_spectro[1:3]*0.1,zans.z, $
      tmpkcorrect,coeffs=coeffs,/stddev, $
      filterlist=['sdss_g0.par','sdss_r0.par','sdss_i0.par'],rmatrix=rmatrix, $
      zvals=zvals, band_shift=0.1
    kcorrect[*,ifiberid]=tmpkcorrect
    lambda=fltarr(n_elements(wave)+1L)
    dwave=alog(wave[1])-alog(wave[0])
    lambda[0:n_elements(wave)-1L]=exp(alog(wave)-0.5*dwave)
    lambda[n_elements(wave)]=exp(alog(wave[n_elements(wave)-1L])+0.5*dwave)
    lambda=lambda/(1.+zans.z)
    vmatrix=flux
    kcorrect,dum1,dum2,zans.z, coeffs=[1.],vmatrix=vmatrix,lambda=lambda, $
      tmpkcorrect,filterlist=['sdss_g0.par','sdss_r0.par','sdss_i0.par'], $
      band_shift=0.1
    kcorrect_real[*,ifiberid]=tmpkcorrect
    print,zans.z
    print,minmax(wave)
    print,kcorrect[*,ifiberid]
    print,kcorrect_real[*,ifiberid]
endfor

stop

end
