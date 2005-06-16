pro test_spec, vname=vname

postcat=hogg_mrdfits(vagc_name('post_catalog', sample='sample15', $
                               letter='bsafe', post='1'), 1, nrow=28800)
postcat=postcat[where(postcat.z gt 0.03 and $
                      postcat.z lt 0.10 and $
                      postcat.m lt 16.4)]
num=100
ii=shuffle_indx(n_elements(postcat), num_sub=num)
postcat=postcat[ii]
sp=mrdfits(vagc_name('object_sdss_spectro'),1,row=postcat.object_position, $
          columns=['fiberid', 'plate', 'mjd'])

hd=fltarr(num)
hderr=fltarr(num)
d4000=fltarr(num)
d4000err=fltarr(num)
hdm=fltarr(num)
d4000m=fltarr(num)

for i=0L, num-1L do begin
    range=[3700., 4300.]
    fit_sdss_fiber, sp[i].plate, sp[i].fiberid, $
      mjd=sp[i].mjd, loglam=loglam, flux=flux, model=model, ivar=ivar, $
      vname=vname, /usev, cmodel=cmodel, range=range
    d4000[i]=get_d4000n(loglam, flux, ivar, ivar=d4000ivar)
    d4000err[i]=1./sqrt(d4000ivar)
    hdstruct=get_line(loglam, flux,  ivar,lname='HDELTA')
    if(n_tags(hdstruct) gt 0) then begin
        hd[i]=hdstruct.hdelta_eqw
        hderr[i]=1./sqrt(hdstruct.hdelta_eqw_ivar)
    endif
    d4000m[i]=get_d4000n(loglam, model)
    hdmstruct=get_line(loglam, model, lname='HDELTA')
    if(n_tags(hdmstruct) gt 0) then begin
        hdm[i]=hdmstruct.hdelta_eqw
    endif
endfor

stop

end
