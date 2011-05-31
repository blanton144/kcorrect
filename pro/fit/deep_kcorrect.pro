;+
; NAME:
;   deep_kcorrect
; PURPOSE:
;   calculate K-corrections for standard DEEP input (from BRI to NUB)
; CALLING SEQUENCE:
;   kcorrect= deep_kcorrect(redshift [, nmgy=, ivar=, mag=, err=, $
;                           zcat=, /sdss, band_shift=, chi2=, rmaggies=, $
;                           omaggies=, vname=, oivar=, mass=, mtol=, $
;                           absmag=, amivar=, omega0=, omegal0= ])
; INPUTS:
;   redshift - [N] redshifts
;   zcat - [N] DEEP zcat-style structure, containing:
;                  .MAGB
;                  .MAGR
;                  .MAGI
;                  .SFD_EBV
;   nmgy, ivar - [3, N] nanomaggies, Galactic-reddening corrected, and inverse
;                variance of same
;   mag, err - [3, N] Pogson magnitudes, Galactic-reddening corrected and
;              errors of same
; OPTIONAL INPUTS:
;   band_shift    - blueshift of output bandpasses (to get ^{z}b
;                   type bands) [default 0.]
;   vname - name of fit to use (defaults to 'default')
;   omega0, omegal0 - cosmological parameters for calculating distance
;                     moduli [default 0.3, 0.7]
;   filterlist - list of filters to K-correct to
; OPTIONAL KEYWORDS:
;   /closest - use closest bands for K-corrections
;   /silent - shut up
; OUTPUTS:
;   kcorrect - [3, N] K-corrections from BRI to NUV, U, and B (or NUV,
;              u, g if /sdss is set) satisfying
;                m_R = M_Q + DM(z) + K_QR(z)
;              based on the best fit sum of templates. All magnitudes
;              are AB.
;   mtol - [3, N] mass-to-light ratios from model in each output band
;   mass - [N] total current stellar mass from model in each band
;   intsfh - [N] total integrated star formation history
;   absmag - [3, N] absolute magnitude (for missing data, substitutes
;            model fit) in each output band
;   amivar - [3, N] inverse variance of absolute magnitude (for
;            missing data = 0) in each output band
; OPTIONAL OUTPUTS:
;   coeffs - coefficients of fit
;   chi2 - chi^2 of fit
;   rmaggies - [3, N] reconstructed maggies from the fit (BRI)
;   omaggies, oivar - [3, N] maggies and inverse variances used for fit
;                           (after extinction, AB correction, etc)
;                           (BRI)
;   obands - [3, N] which input bands the K-corrections refer to 
; COMMENTS:
;   This is a simple wrapper on kcorrect.pro.  It keeps a version of
;   rmatrix and zvals in memory to save time, recalculating them each
;   time you change band_shift.
;
;   You must specify nmgy,ivar OR mag,err OR zcat. If
;   nmgy or mag, make sure they are AB calibrated and Galactic
;   extinction corrected.
;
;   Uses deep_to_maggies to convert zcat structure to Galactic
;   extinction corrected maggies with errors.
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
function deep_kcorrect, redshift, nmgy=nmgy, ivar=ivar, mag=mag, err=err, $
                        zcat=zcat, band_shift=in_band_shift, chi2=chi2, $
                        coeffs=coeffs, rmaggies=rmaggies, omaggies=omaggies, $
                        oivar=oivar, vname=in_vname, mass=mass, mtol=mtol, $
                        absmag=absmag, amivar=amivar, sdss=sdss, $
                        rmatrix=rmatrix, closest=closest, obands=obands, $
                        omega0=omega0, omegal0=omegal0, $
                        filterlist=in_filterlist, silent=silent, intsfh=intsfh

common com_deep_kcorrect, out_rmatrix, out_zvals, band_shift, $
  deep_rmatrix, deep_zvals, vname, out_filterlist

if(n_params() lt 1 OR $
   (((keyword_set(nmgy) eq 0 OR keyword_set(ivar) eq 0)) AND $
    ((keyword_set(mag) eq 0 OR keyword_set(err) eq 0)) AND $
    (n_tags(zcat) eq 0) AND $
    (n_elements(coeffs) eq 0))) $
  then begin
    doc_library, 'deep_kcorrect'
    return, -1
endif 

new_out_filterlist=['galex_NUV.par', $
                    'bessell_U.par', $
                    'bessell_B.par']
if(keyword_set(sdss)) then $
  new_out_filterlist=['galex_NUV.par', $
                      'sdss_u0.par', $
                      'sdss_g0.par']
if(keyword_set(in_filterlist)) then $
  new_out_filterlist=in_filterlist

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
        deep_rmatrix=0
        deep_zvals=0
        out_rmatrix=0
        out_zvals=0
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
    if(n_elements(out_filterlist) ne n_elements(new_out_filterlist)) then begin
        out_rmatrix=0
        out_zvals=0
    endif
    for i=0L, n_elements(out_filterlist)-1L do begin
        if(new_out_filterlist[i] ne out_filterlist[i]) then begin
            out_rmatrix=0
            out_zvals=0
        endif
    endfor
endif 
out_filterlist=new_out_filterlist

if(keyword_set(mag) or keyword_set(nmgy)) then begin
    mgy=fltarr(3, n_elements(redshift))
    mgy_ivar=fltarr(3, n_elements(redshift))
    if(n_elements(mag) gt 0) then begin
        mgy[*,*]=10.^(-0.4*mag)
        mgy_ivar[*,*]=1./(0.4*alog(10.)*mgy*err)^2.
    endif
    if(n_elements(nmgy) gt 0) then begin
        mgy[*,*]=nmgy*1.e-9
        mgy_ivar[*,*]=ivar*1.e+18
    endif
endif

if(n_tags(zcat) gt 0) then $
  deep_to_maggies, zcat, mgy, mgy_ivar

;; call kcorrect
deep_filterlist=['deep_B.par', 'deep_R.par', 'deep_I.par']
kcorrect, mgy, mgy_ivar, redshift, kcdum, band_shift=band_shift, $
  rmatrix=deep_rmatrix, zvals=deep_zvals, coeffs=coeffs, rmaggies=rmaggies, $
  vname=vname, mass=mass, mtol=mtol, absmag=absmag, amivar=amivar, $
  filterlist=deep_filterlist, silent=silent, intsfh=intsfh

; calculate the preliminaries
if(NOT keyword_set(out_rmatrix) OR NOT keyword_set(out_zvals)) then begin
    if(NOT keyword_set(vmatrix) OR NOT keyword_set(lambda)) then $
      k_load_vmatrix, vmatrix, lambda, vfile=vfile, lfile=lfile, $
      vpath=vpath, vname=vname
    k_projection_table,out_rmatrix,vmatrix,lambda,out_zvals,out_filterlist, $ 
      zmin=zmin,zmax=zmax,nz=nz,filterpath=filterpath, silent=silent
endif

k_reconstruct_maggies,coeffs,replicate(band_shift,n_elements(redshift)), $
  reconstruct_maggies,rmatrix=out_rmatrix,zvals=out_zvals
reconstruct_maggies=reconstruct_maggies/(1.+band_shift)

obands=lindgen(n_elements(out_filterlist))#replicate(1L, n_elements(redshift))

if(keyword_set(closest)) then begin
    lambda_in=k_lambda_eff(filterlist=deep_filterlist)
    lambda_out=k_lambda_eff(filterlist=out_filterlist, band_shift=band_shift)
    for i=0L, n_elements(redshift)-1L do begin
        for j=0L, n_elements(lambda_out)-1L do begin
            dmin=min(abs(lambda_in/(1.+redshift[i])-lambda_out[j]), imin)
            obands[j, i]= imin
        endfor
    endfor
endif

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

rmatrix=out_rmatrix

return, kcorrect

end
