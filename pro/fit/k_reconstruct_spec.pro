;+
; NAME:
;   k_reconstruct_spec
; PURPOSE:
;   Reconstruct galaxy rest-frame spectrum given a fit
; CALLING SEQUENCE:
;   k_reconstruct_spec, coeffs, loglam, flux [, vname=, vdisp=, $
;        /nolines, /noextinct, /init, nt= ]
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
; COMMENTS:
;   If vdisp, /nolines, and /noextinct 
;   If lines are included, they are always smoothed at 300 km/s vdisp
;   Bases fit on file:
;      $KCORRECT_DIR/data/templates/k_nmf_tspec.[vname].fits
;   If such a file doesn't exist, builds it from:
;      $KCORRECT_DIR/data/templates/k_nmf_mmatrix.[vname].fits
;      $KCORRECT_DIR/data/templates/k_nmf_rawspec.[vname].fits
;      $KCORRECT_DIR/data/templates/k_nmf_soln.[vname].fits
; REVISION HISTORY:
;   21-Apr-2005  Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_reconstruct_spec,coeffs, ologlam, flux, vname=in_vname, vdisp=vdisp, $
                       nolines=nolines, noextinct=noextinct, nt=out_nt, $
                       init=init, reset=reset

common com_krs, vname, nt, loglam, lspec_v300, $
  tspec_v300, tspec_v300_nl, tspec_v300_nl_nd, tspec_v300_nd, $
  tspec_v0,   tspec_v0_nl,   tspec_v0_nl_nd,   tspec_v0_nd

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
    tspecfile=getenv('KCORRECT_DIR')+'/data/templates/k_nmf_spec.'+ $
      vname+'.fits'
    if(file_test(tspecfile) eq 0 OR keyword_set(reset) eq 1) then begin
        spec=mrdfits(getenv('KCORRECT_DIR')+'/data/templates/k_nmf_mmatrix.'+ $
                     vname+'.fits', 0, hdr)
        nel=long(sxpar(hdr,'NEL'))
        nspec=long(sxpar(hdr,'NSPEC'))
        lambda=mrdfits(getenv('KCORRECT_DIR')+ $
                       '/data/templates/k_nmf_mmatrix.'+ $
                       vname+'.fits', 1)
        dust=mrdfits(getenv('KCORRECT_DIR')+'/data/templates/k_nmf_mmatrix.'+ $
                     vname+'.fits', 2)
        rawspec=mrdfits(getenv('KCORRECT_DIR')+ $
                        '/data/templates/k_nmf_rawspec.'+ $
                        vname+'.fits', 0)
        templates=mrdfits(getenv('KCORRECT_DIR')+ $
                          '/data/templates/k_nmf_soln.'+ $
                          vname+'.fits', 0)
        
        nl=n_elements(lambda)
        nb=n_elements(spec)/nl
        ns=n_elements(dust)
        nt=n_elements(templates)/nb
        
        loglam=alog10(lambda[0:nspec-1L])
        absrc=3.631*2.99792*1.e-2/lambda^2
        for i=0L, nb-1L do spec[*,i]=spec[*,i]*absrc
        
        tspec_v300=spec[0L:nspec-1L,*]#templates
        tspec_v300_nl= $
          spec[0:nspec-1L,0:n_elements(dust)-1L]# $
          templates[0:n_elements(dust)-1L,*]
        tspec_v0_nl=rawspec[0L:nspec-1L,*]#templates[0:n_elements(dust)-1L,*]
        
        ;; assumes dust extinction = 0 is in first section
        ii=where(dust.tauv eq 0, nii)
        nd=ns/nii
        templates_nd= $
          total(reform(templates[0:n_elements(dust)-1L,*], nii, nd,nt), 2)
        tspec_v300_nl_nd=spec[0:nspec-1L,0:nii-1L]#templates_nd
        tspec_v0_nl_nd=rawspec[0:nspec-1L,0:nii-1L]#templates_nd
        
        ;; assumes lines are right after spectrum
        lspec_v300=spec[0L:nspec-1L,n_elements(dust): $
                        n_elements(dust)+nel-1L]# $
          templates[n_elements(dust):n_elements(dust)+nel-1L, *]
        tspec_v0=tspec_v0_nl+lspec_v300
        tspec_v300_nd=tspec_v300_nl_nd+lspec_v300
        tspec_v0_nd=tspec_v0_nl_nd+lspec_v300

        hdr=['']
        sxaddpar, hdr, 'NT', nt, 'number of templates'
        mwrfits, loglam, tspecfile, hdr, /create
        mwrfits, tspec_v0, tspecfile
        mwrfits, tspec_v0_nl, tspecfile
        mwrfits, tspec_v0_nd, tspecfile
        mwrfits, tspec_v0_nl_nd, tspecfile
        mwrfits, tspec_v300, tspecfile
        mwrfits, tspec_v300_nl, tspecfile
        mwrfits, tspec_v300_nd, tspecfile
        mwrfits, tspec_v300_nl_nd, tspecfile
        mwrfits, lspec_v300, tspecfile
    endif else begin
        loglam=mrdfits(tspecfile, 0, hdr)
        nt=long(sxpar(hdr,'NT'))
        tspec_v0=mrdfits(tspecfile, 1)
        tspec_v0_nl=mrdfits(tspecfile, 2)
        tspec_v0_nd=mrdfits(tspecfile, 3)
        tspec_v0_nd_nl=mrdfits(tspecfile, 4)
        tspec_v300=mrdfits(tspecfile, 5)
        tspec_v300_nl=mrdfits(tspecfile, 6)
        tspec_v300_nd=mrdfits(tspecfile, 7)
        tspec_v300_nd_nl=mrdfits(tspecfile, 8)
        lspec_v300=mrdfits(tspecfile, 9)
    endelse
endif

out_nt=nt
ologlam=loglam

if(keyword_set(init)) then return

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
            sigma=sqrt(vdisp)/(2.99792e+5*alog(10.))
            flux=gauss_smooth(loglam,flux,sigma,loglam)
        endif
        return
    endif
    
    if(keyword_set(nolines) eq 1 and $
       keyword_set(noextinct) eq 0) then begin
        flux=tspec_v0_nl#coeffs 
        if(vdisp gt 0.) then begin
            sigma=sqrt(vdisp)/(2.99792e+5*alog(10.))
            flux=gauss_smooth(loglam,flux,sigma,loglam)
        endif
        return
    endif
    
    if(keyword_set(nolines) eq 0 and $
       keyword_set(noextinct) eq 0) then begin
        flux=tspec_v0_nl#coeffs 
        if(vdisp gt 0.) then begin
            sigma=vdisp/(2.99792e+5*alog(10.))
            flux=gauss_smooth(loglam,flux,sigma,loglam)
        endif
        flux=flux+lspec_v300#coeffs
        return
    endif
    
    if(keyword_set(nolines) eq 0 and $
       keyword_set(noextinct) eq 1) then begin
        flux=tspec_v0_nl_nd#coeffs 
        if(vdisp gt 0.) then begin
            sigma=sqrt(vdisp)/(2.99792e+5*alog(10.))
            flux=gauss_smooth(loglam,flux,sigma,loglam)
        endif
        flux=flux+lspec_v300#coeffs
        return
    endif
endelse

end
;------------------------------------------------------------------------------
