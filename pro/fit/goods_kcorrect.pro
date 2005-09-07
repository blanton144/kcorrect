;+
; NAME:
;   goods_kcorrect
; PURPOSE:
;   calculate K-corrections for GOODS catalog (BVizJHK)
; CALLING SEQUENCE:
;   kcorrect= goods_kcorrect(redshift [, nmgy=, ivar=, mag=, err=, $
;                            goods=, /sdss, band_shift=, chi2=, rmaggies=, $
;                            omaggies=, vname=, oivar=, mass=, mtol=, $
;                            absmag=, amivar=, omega0=, omegal0= ])
; INPUTS:
;   redshift - [N] redshifts
;   goods - [N] GOODS-style structure, containing:
;               .RA (J2000 degrees)
;               .DEC (J2000 degrees)
;               .BMAG_MAGAUTO
;               .BMAGERR_MAGAUTO
;               .VMAG_MAGAUTO
;               .VMAGERR_MAGAUTO
;               .IMAG_MAGAUTO
;               .IMAGERR_MAGAUTO
;               .ZMAG_MAGAUTO
;               .ZMAGERR_MAGAUTO
;               .JMAG_MAGAUTO
;               .JMAGERR_MAGAUTO
;               .HMAG_MAGAUTO
;               .HMAGERR_MAGAUTO
;               .KMAG_MAGAUTO
;               .KMAGERR_MAGAUTO
;   nmgy, ivar - [7, N] nanomaggies, Galactic-reddening corrected, and inverse
;                variance of same
;   mag, err - [7, N] Pogson magnitudes, Galactic-reddening corrected and
;              errors of same
; OPTIONAL INPUTS:
;   band_shift    - blueshift of output bandpasses (to get ^{z}b
;                   type bands) [default 0.]
;   vname - name of fit to use [defaults to 'default']
;   omega0, omegal0 - cosmological parameters for calculating distance
;                     moduli [default 0.3, 0.7]
;   filterlist - [Nf] output filter list; default is:
;                   ['galex_NUV.par', $
;                    'bessell_U.par', $
;                    'bessell_B.par', $
;                    'bessell_V.par', $
;                    'bessell_R.par', $
;                    'bessell_I.par']
;                unless /sdss is set, in which case: 	
;                   ['galex_NUV.par', $
;                    'sdss_u0.par', $
;                    'sdss_g0.par', $
;                    'sdss_r0.par', $
;                    'sdss_i0.par', $
;                    'sdss_z0.par']
; OUTPUTS:
;   kcorrect - [Nf, N] K-corrections from BVizJHK to NUBVRI (or Nugriz 
;              /sdss is set). The closest input band is used for each
;              output band (which ones are decided upon is output in
;              "obands"). K-corrections satisfy
;                m_R = M_Q + DM(z) + K_QR(z)
;              based on the best fit sum of templates. All magnitudes
;              are AB. 
;   mtol - [Nf, N] mass-to-light ratios from model in each output band
;   mass - [N] total mass from model in each band
;   absmag - [Nf, N] absolute magnitude (for missing data, substitutes
;            model fit) in each output band
;   amivar - [Nf, N] inverse variance of absolute magnitude (for
;            missing data = 0) in each output band
; OPTIONAL OUTPUTS:
;   coeffs - coefficients of fit
;   chi2 - chi^2 of fit
;   rmaggies - [7, N] reconstructed maggies from the fit (BVizJHK)
;   omaggies, oivar - [7, N] maggies and inverse variances used for fit
;                           (after extinction correction, etc)
;                           (BVizJHK)
;   obands - [Nf, N] which input bands the K-corrections refer to 
;   appm - [Nf, N] apparent magnitudes in output bands
; COMMENTS:
;   This is a simple wrapper on kcorrect.pro.  It keeps a version of
;   rmatrix and zvals in memory to save time, recalculating them each
;   time you change band_shift.
;
;   You must specify nmgy,ivar OR mag,err OR goods. If
;   nmgy or mag, make sure they are AB calibrated and Galactic
;   extinction corrected.
;
;   Uses goods_to_maggies to convert goods structure to Galactic
;   extinction corrected maggies with errors. This ignores the H band,
;   which we think is crap.
;
;   For v4_0b templates and later, coefficients are in units of: 
; 
;     1 solar mass / (D/10pc)^2 
;
;   That is, sum the coefficients and multiply by (D/10pc)^2 to get
;   TOTAL INTEGRATED STAR FORMATION. (In fact, for Omega0=0.3 and
;   OmegaL0=0.7, this is what the "mass" keyword returns). Note that
;   the total integrated star formation DIFFERS from the current
;   stellar mass --- which is returned in the mass and mtol variables.
; REVISION HISTORY:
;   07-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
function goods_kcorrect, redshift, nmgy=nmgy, ivar=ivar, mag=mag, err=err, $
                         goods=goods, band_shift=in_band_shift, chi2=chi2, $
                         coeffs=coeffs, rmaggies=rmaggies, omaggies=omaggies, $
                         oivar=oivar, vname=in_vname, mass=mass, mtol=mtol, $
                         absmag=absmag, amivar=amivar, sdss=sdss, $
                         rmatrix=rmatrix, obands=obands, omega0=omega0, $
                         omegal0=omegal0, filterlist=filterlist, appm=appm, $
                         useh=useh

common com_goods_kcorrect, out_rmatrix, out_zvals, band_shift, $
  goods_rmatrix, goods_zvals, vname, out_filterlist

if(n_params() lt 1 OR $
   (((keyword_set(nmgy) eq 0 OR keyword_set(ivar) eq 0)) AND $
    ((keyword_set(mag) eq 0 OR keyword_set(err) eq 0)) AND $
    (n_tags(goods) eq 0))) $
  then begin
    doc_library, 'goods_kcorrect'
    return, -1
endif 

new_out_filterlist=['galex_NUV.par', $
                    'bessell_U.par', $
                    'bessell_B.par', $
                    'bessell_V.par', $
                    'bessell_R.par', $
                    'bessell_I.par']
if(keyword_set(sdss)) then $
  new_out_filterlist=['galex_NUV.par', $
                      'sdss_u0.par', $
                      'sdss_g0.par', $
                      'sdss_r0.par', $
                      'sdss_i0.par', $
                      'sdss_z0.par']
if(keyword_set(filterlist)) then $
  new_out_filterlist=filterlist

if(n_elements(in_vname) gt 0) then begin
    use_vname=in_vname
endif else begin
    if(keyword_set(lrg)) then $
      use_vname='lrg1' $
    else $
      use_vname='default'
endelse
if(n_elements(vname) gt 0) then begin
    if(vname ne use_vname) then begin
        rmatrix=0
        zvals=0
        goods_rmatrix=0
        goods_zvals=0
        out_rmatrix=0
        out_zvals=0
        vmatrix=0
        lambda=0
    endif
endif
vname=use_vname

;; need to reset rmatrix if band_shift changes
if(n_elements(in_band_shift) gt 0) then $
  use_band_shift=in_band_shift $
else $
  use_band_shift=0. 
if(n_elements(band_shift) gt 0) then begin
    if(band_shift ne use_band_shift) then begin
        out_rmatrix=0
        out_zvals=0
    endif
endif
band_shift=use_band_shift

;; need to reset rmatrix if filterlist changes
if(n_elements(out_filterlist) gt 0) then begin
    if(n_elements(out_filterlist) eq n_elements(new_out_filterlist)) then begin
        for i=0L, n_elements(out_filterlist)-1L do begin
            if(new_out_filterlist[i] ne out_filterlist[i]) then begin
                out_rmatrix=0
                out_zvals=0
            endif
        endfor
    endif else begin
        out_rmatrix=0
        out_zvals=0
    endelse
endif 
out_filterlist=new_out_filterlist

mgy=fltarr(7, n_elements(redshift))
mgy_ivar=fltarr(7, n_elements(redshift))
if(n_elements(mag) gt 0) then begin
    mgy[*,*]=10.^(-0.4*mag)
    mgy_ivar[*,*]=1./(0.4*alog(10.)*mgy*err)^2.
endif
if(n_elements(nmgy) gt 0) then begin
    mgy[*,*]=nmgy*1.e-9
    mgy_ivar[*,*]=ivar*1.e+18
endif

if(n_tags(goods) gt 0) then $
  goods_to_maggies, goods, mgy, mgy_ivar, useh=useh

;; call kcorrect
goods_filterlist=['goods_acs_f435w.par', $
                  'goods_acs_f606w.par', $
                  'goods_acs_f775w.par', $
                  'goods_acs_f850lp.par', $
                  'goods_J_isaac_etc.par', $
                  'goods_H_isaac_etc.par', $
                  'goods_Ks_isaac_etc.par']
kcorrect, mgy, mgy_ivar, redshift, kcdum, band_shift=band_shift, $
  rmatrix=goods_rmatrix, zvals=goods_zvals, coeffs=coeffs, rmaggies=rmaggies, $
  vname=vname, mass=mass, mtol=mtol, absmag=absmag, amivar=amivar, $
  filterlist=goods_filterlist

; calculate the preliminaries
if(NOT keyword_set(out_rmatrix) OR NOT keyword_set(out_zvals)) then begin
    if(NOT keyword_set(vmatrix) OR NOT keyword_set(lambda)) then $
      k_load_vmatrix, vmatrix, lambda, vfile=vfile, lfile=lfile, $
      vpath=vpath, vname=vname
    k_projection_table,out_rmatrix,vmatrix,lambda,out_zvals,out_filterlist, $ 
      zmin=zmin,zmax=zmax,nz=nz,filterpath=filterpath
endif

k_reconstruct_maggies,coeffs,replicate(band_shift,n_elements(redshift)), $
  reconstruct_maggies,rmatrix=out_rmatrix,zvals=out_zvals
reconstruct_maggies=reconstruct_maggies/(1.+band_shift)

obands=lindgen(n_elements(out_filterlist))#replicate(1L, n_elements(redshift))

lambda_in=k_lambda_eff(filterlist=goods_filterlist)
lambda_out=k_lambda_eff(filterlist=out_filterlist, band_shift=band_shift)
for i=0L, n_elements(redshift)-1L do begin
    for j=0L, n_elements(lambda_out)-1L do begin
        dmin=min(abs(lambda_in/(1.+redshift[i])-lambda_out[j]), imin)
        obands[j, i]= imin
    endfor
endfor

kcorrect=fltarr(n_elements(out_filterlist), n_elements(redshift))
for i=0L, n_elements(redshift)-1L do $
  for j=0L, n_elements(out_filterlist)-1L do $
  kcorrect[j,i]=reconstruct_maggies[j,i]/ rmaggies[obands[j,i],i]
kcorrect=2.5*alog10(kcorrect)

if(arg_present(omaggies) or arg_present(absmag) gt 0) then $
  omaggies=mgy
if(arg_present(oivar) gt 0 or arg_present(absmag) gt 0) then $
  oivar=mgy_ivar

if(arg_present(absmag)) then begin
    absmag=fltarr(n_elements(out_filterlist), n_elements(redshift))
    amivar=fltarr(n_elements(out_filterlist), n_elements(redshift))
    for i=0L, n_elements(out_filterlist)-1L do $
      absmag[i,*]=-2.5*alog10(reconstruct_maggies[i,*])- $
      lf_distmod(redshift, omega0=omega0, omegal0=omegal0)
    for j=0L, n_elements(redshift)-1L do begin
        igood=where(oivar[obands[*,j],j] gt 0. and $
                    omaggies[obands[*,j],j] gt 0., ngood)
        if(ngood gt 0) then begin
            absmag[igood,j]=-2.5*alog10(omaggies[obands[igood,j],j])- $
              lf_distmod(redshift[j],omega0=omega0,omegal0=omegal0)- $
              kcorrect[igood,j]
            amivar[igood,j]=omaggies[obands[igood,j],j]^2* $
              oivar[obands[igood,j],j]*(0.4*alog(10.))^2
        endif
    endfor
endif

if(arg_present(appm)) then begin
    ;; get apparent magnitudes in output bands 
    k_reconstruct_maggies,coeffs,redshift, $
      app_maggies,rmatrix=out_rmatrix,zvals=out_zvals

    ;; get comparisons
    appm=fltarr(n_elements(out_filterlist), n_elements(redshift))
    appm_obands=lindgen(n_elements(out_filterlist))
    lambda_in=k_lambda_eff(filterlist=goods_filterlist)
    lambda_out=k_lambda_eff(filterlist=out_filterlist, band_shift=band_shift)
    for j=0L, n_elements(lambda_out)-1L do begin
        dmin=min(abs(lambda_in-lambda_out[j]), imin)
        appm_obands[j]= imin
    endfor

    appmcorrect=fltarr(n_elements(out_filterlist), n_elements(redshift))
    for j=0L, n_elements(out_filterlist)-1L do $
      appmcorrect[j,*]=app_maggies[j,*]/rmaggies[appm_obands[j],*]
    appmcorrect=2.5*alog10(appmcorrect)

    appm=fltarr(n_elements(out_filterlist), n_elements(redshift))
    for i=0L, n_elements(out_filterlist)-1L do $
      appm[i,*]=-2.5*alog10(app_maggies[i,*])
    for j=0L, n_elements(redshift)-1L do begin
        igood=where(oivar[appm_obands,j] gt 0. and $
                    omaggies[appm_obands,j] gt 0., ngood)
        if(ngood gt 0) then $
          appm[igood,j]=-2.5*alog10(omaggies[appm_obands[igood],j])- $
          appmcorrect[igood,j]
    endfor
    
endif

rmatrix=out_rmatrix

return, kcorrect

end
