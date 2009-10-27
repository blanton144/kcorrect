;+
; NAME:
;   k_reconstruct_spec
; PURPOSE:
;   Reconstruct galaxy rest-frame spectrum given a fit
; CALLING SEQUENCE:
;   k_reconstruct_spec, coeffs, loglam, flux [, vname=, vdisp=, $
;        /nolines, /noextinct, /init, nt=, mass=, metallicity=, $
;        b300=, b1000=]
; INPUTS:
;   coeffs - [NT, NGALS] coefficients 
; OUTPUTS:
;   loglam - [NL, NGALS] wavelengths
;   flux - [NL, NGALS] fluxes (ergs cm-2 s-1 A-1)
; OPTIONAL KEYWORDS
;   /nolines - do not include lines
;   /noextinct - show unextincted spectra
;   /init - don't make a spectrum, just initialize
;   /reset - create the k_nmf_tspec.fits file if it doesn't exist
; OPTIONAL INPUTS:
;   vname - name of fit to use (default 'default')
;   vdisp - smooth with this velocity dispersion
; OPTIONAL OUTPUTS:
;   nt - total number of templates
;   mass, metallicity - properties of template fit;
;                       mass is current stellar mass and is in units of 
;                       1 solar mass / (D/10pc)^2
;   b300 - star-formation within last 300Myrs relative to average
;          star-formation rate
;   b1000 - star-formation within last 1Gyrs relative to average
;           star-formation rate
; COMMENTS:
;   If coeffs are standard, returns units of erg/cm^2/s/A
;   If vdisp, /nolines, and /noextinct 
;   If lines are included, they are always smoothed at 300 km/s vdisp
;   Bases fit on file:
;      $KCORRECT_DIR/data/templates/k_nmf_derived.[vname].fits
; REVISION HISTORY:
;   21-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_reconstruct_spec,coeffs, ologlam, flux, vname=in_vname, vdisp=vdisp, $
                       nolines=nolines, noextinct=noextinct, nt=out_nt, $
                       init=init, reset=reset, mass=mass, $
                       metallicity=metallicity, b300=b300, b1000=b1000

common com_krs, vname, nt, loglam, lspec_v300, tmass, tmetallicity, $
  tspec_v300, tspec_v300_nl, tspec_v300_nl_nd, tspec_v300_nd, $
  tspec_v0,   tspec_v0_nl,   tspec_v0_nl_nd,   tspec_v0_nd, tmass300, $
  tmremain, tmass1000

; Need at least 3 parameters
if (N_params() LT 1) then begin
    doc_library, 'k_reconstruct_spec'
    return
endif

if (NOT keyword_set(in_vname)) then $
  in_vname='default'

newname=1
if(n_elements(vname) gt 0) then $
  if(vname eq in_vname) then newname=0
vname=in_vname

if(newname) then begin
    tspecfile= $
      getenv('KCORRECT_DIR')+'/data/templates/k_nmf_derived.'+vname+'.fits'
    hdr=headfits(tspecfile)
    nt=long(sxpar(hdr,'NT'))
    loglam=alog10(mrdfits(tspecfile, 11))
    tspec_v0=mrdfits(tspecfile, 1)
    tspec_v0_nl=mrdfits(tspecfile, 2)
    tspec_v0_nd=mrdfits(tspecfile, 3)
    tspec_v0_nd_nl=mrdfits(tspecfile, 4)
    tspec_v300=mrdfits(tspecfile, 5)
    tspec_v300_nl=mrdfits(tspecfile, 6)
    tspec_v300_nd=mrdfits(tspecfile, 7)
    tspec_v300_nd_nl=mrdfits(tspecfile, 8)
    lspec_v300=mrdfits(tspecfile, 9)
    tmass=mrdfits(tspecfile, 17)
    tmetallicity=mrdfits(tspecfile, 18)
    tmass300=mrdfits(tspecfile, 19)
    tmass1000=mrdfits(tspecfile, 20)
    tmremain=mrdfits(tspecfile, 24)
endif

out_nt=nt
ologlam=loglam

if(keyword_set(init)) then return

mass=total(tmremain*coeffs)
b300=total(tmass300*coeffs)/total(tmass*coeffs)
b1000=total(tmass1000*coeffs)/total(tmass*coeffs)
metallicity=total(tmremain*tmetallicity*coeffs)/mass

if(n_elements(vdisp) eq 0) then begin

    if(keyword_set(nolines) eq 1 and $
       keyword_set(noextinct) eq 1) then begin
        flux=tspec_v300_nl_nd#coeffs 
        return
    endif
      
    if(keyword_set(nolines) eq 1 and $
       keyword_set(noextinct) eq 0) then begin
        flux=tspec_v300_nl#coeffs 
        return
    endif
      
    if(keyword_set(nolines) eq 0 and $
       keyword_set(noextinct) eq 0) then begin
        flux=tspec_v300#coeffs 
        return
    endif
      
    if(keyword_set(nolines) eq 0 and $
       keyword_set(noextinct) eq 1) then begin
        flux=tspec_v300_nd#coeffs 
        return
    endif

endif else begin

    if(keyword_set(nolines) eq 1 and $
       keyword_set(noextinct) eq 1) then begin
        flux=tspec_v0_nl_nd#coeffs 
        if(vdisp gt 0.) then begin
            flux=k_smooth(loglam,flux,(vdisp))
        endif
        return
    endif
    
    if(keyword_set(nolines) eq 1 and $
       keyword_set(noextinct) eq 0) then begin
        flux=tspec_v0_nl#coeffs 
        if(vdisp gt 0.) then begin
            flux=k_smooth(loglam,flux,(vdisp))
        endif
        return
    endif
    
    if(keyword_set(nolines) eq 0 and $
       keyword_set(noextinct) eq 0) then begin
        flux=tspec_v0_nl#coeffs 
        if(vdisp gt 0.) then begin
            flux=k_smooth(loglam,flux,(vdisp))
        endif
        flux=flux+lspec_v300#coeffs
        return
    endif
    
    if(keyword_set(nolines) eq 0 and $
       keyword_set(noextinct) eq 1) then begin
        flux=tspec_v0_nl_nd#coeffs 
        if(vdisp gt 0.) then begin
            flux=k_smooth(loglam,flux,(vdisp))
        endif
        flux=flux+lspec_v300#coeffs
        return
    endif
endelse

end
;------------------------------------------------------------------------------
