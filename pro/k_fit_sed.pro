;+
; NAME:
;   k_fit_sed
;
; PURPOSE:
;   Fit SED using k_fit_coeffs and k_fit_templates iteratively
;
; CALLING SEQUENCE:
;   k_fit_sed, galaxy_flux, galaxy_invvar, galaxy_z, templatelist, $
;     filterlist, coeff, ematrix, bmatrix, bflux, lambda, $
;     [smoothtemplate=, sublmin=, sublmax=, vmatrix=, $
;     preset_ematrix=, maxiter=, nt=, model_flux=, plotfluxes=]
;
; INPUTS:
;   galaxy_flux   - flux in each band for each galaxy [N_band, N_gal]
;   galaxy_invvar - errors in each band for each galaxy [N_band, N_gal]
;   galaxy_z      - redshift for each galaxy [N_gal]
;   templatelist  - list of templates specifying full SED space
;   filterlist    - list of files with filter information [N_band]
;
; OPTIONAL INPUTS:
;   smoothtemplate   - smooth templates by this many angstroms
;   sublmin, sublmax - range for estimate of flux and orthogonalization
;                      in Angstroms [3500, 7500]
;   preset_ematrix   - initial set of eigentemplates to use [N_dim, N_t]
;                      (defaults to unity on the diagonal)
;   maxiter          - maximum number of iterations (default 10)
;   nt               - number of eigentemplates (default N_band-1)
;   model_flux       - reconstructed fluxes
;   plot_fluxes      - set to plot pluxes as fit is done for debugging
;
; OUTPUTS:
;   coeff    - output coefficients [N_template, N_gal]
;   ematrix       - eigentemplates [N_dim, N_template]
;   bmatrix       - orthogonal templates spanning SED space [N_lambda, N_dim]
;   bflux         - flux in the sublmin,sublmax range for each template [N_dim]
;   lambda        - wavelengths for orthogonal templates [N_lambda]
;
; OPTIONAL INPUT/OUTPUTS:
;   vmatrix          - matrix of initial templates [N_lambda, N_dim]
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;   k_load_templates
;   k_ortho_templates
;   k_create_r
;   k_ortho_etemplates
;   k_fit_coeffs
;   k_fit_templates
;
; REVISION HISTORY:
;   05-Jan-2002  Translated to IDL by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_fit_sed, galaxy_flux, galaxy_invvar, galaxy_z, templatelist, filterlist, coeff, ematrix, bmatrix, bflux, lambda, smoothtemplate=smoothtemplate, sublmin=sublmin, sublmax=sublmax, vmatrix=vmatrix, preset_ematrix=preset_ematrix, maxiter=maxiter, nt=nt, model_flux=model_flux, plotfluxes=plotfluxes, cutlmin=cutlmin, cutlmax=cutlmax, subsmoothtemplate=subsmoothtemplate, subsmoothlimits=subsmoothlimits

; Need at least 10 parameters
if (N_params() LT 10) then begin
    klog, 'Syntax - k_fit_sed, galaxy_flux, galaxy_invvar, galaxy_z, templatelist, $'
    klog, '    filterlist, coeff, ematrix, bmatrix, bflux, lambda, [smoothtemplate=, $'
    klog, '    sublmin=, sublmax=, vmatrix=, preset_ematrix=, maxiter=, nt=, $'
    klog, '    model_flux=, /plotfluxes]'
    return
endif

; Set defaults
if (NOT keyword_set(maxiter)) then maxiter=100l
if (NOT keyword_set(sublmin)) then sublmin=3500.d
if (NOT keyword_set(sublmax)) then sublmax=7500.d
if (NOT keyword_set(cutlmin)) then cutlmin=2000.d
if (NOT keyword_set(cutlmax)) then cutlmax=11500.d

; Find vmatrix and bmatrix
k_load_templates,getenv('KCORRECT_DIR')+'/data/seds/'+templatelist,vmatrix, $
  lambda
nl=long(n_elements(lambda))-1l
nb=long(n_elements(vmatrix))/nl
if(keyword_set(smoothtemplate)) then begin
    sl=(lambda[nl]-lambda[0])/double(nl)
    if(keyword_set(subsmoothtemplate)) then begin 
        ng=long(8*subsmoothtemplate/sl)
        gaussian=exp(-0.5*((dindgen(ng)-0.5*double(ng))/(double(ng)/8))^2)
        gaussian=gaussian/total(gaussian,/double)
        subvmatrix=dblarr(nl,nb)
        for b=0l, nb-1l do begin
            subvmatrix[*,b]=convol(vmatrix[*,b],gaussian,/edge_truncate)
        endfor
    endif
    ng=long(8*smoothtemplate/sl)
    gaussian=exp(-0.5*((dindgen(ng)-0.5*double(ng))/(double(ng)/8))^2)
    gaussian=gaussian/total(gaussian,/double)
    for b=0l, nb-1l do begin
        vmatrix[*,b]=convol(vmatrix[*,b],gaussian,/edge_truncate)
    endfor
    if(keyword_set(subvmatrix)) then begin
        subvmatrix=subvmatrix-vmatrix
        addsub=dblarr(nl,nb)
        for b=0, nb-1 do begin
            indx=where(lambda[0:nl-1] gt subsmoothlimits[0] and $
                       lambda[0:nl-1] lt subsmoothlimits[1])
            addsub[indx,b]=1.d
            indx=where(lambda[0:nl-1] gt subsmoothlimits[0] and $
                       lambda[0:nl-1] lt subsmoothlimits[0]+ $
                       2.*smoothtemplate)
            addsub[indx,b]= $
              0.5*(1.-cos(3.14159*(lambda[indx]-subsmoothlimits[0]) $
                          /(2.*smoothtemplate)))
            indx=where(lambda[0:nl-1] lt subsmoothlimits[1] and $
                       lambda[0:nl-1] gt subsmoothlimits[1]- $
                       2.*smoothtemplate)
            addsub[indx,b]= $
              0.5*(1.-cos(3.14159*(lambda[indx]-subsmoothlimits[1]) $
                          /(2.*smoothtemplate)))
        endfor
        vmatrix=vmatrix+subvmatrix*addsub
    endif
endif
k_ortho_templates,vmatrix,lambda,bmatrix,bflux,sublmin=sublmin,sublmax=sublmax

; Find rmatrix
k_create_r,rmatrix,bmatrix,lambda,zvals,filterlist

; Iterate 
ngalaxy=long(n_elements(galaxy_z))
nk=long(n_elements(galaxy_flux))/ngalaxy
if (NOT keyword_set(nt)) then nt=nk-1l
if (NOT keyword_set(preset_ematrix)) then begin
  ematrix=dblarr(nb,nt)
  indx=lindgen(nt)
	ematrix[indx,indx]=1.d
endif
klog,'start iteration'
for i=0l, maxiter-1l do begin
  k_ortho_etemplates,ematrix,bflux
  k_fit_coeffs,galaxy_flux,galaxy_invvar,galaxy_z,coeff,ematrix=ematrix, $
    zvals=zvals, rmatrix=rmatrix
  k_pca_etemplates,coeff,ematrix,bflux
  k_fit_templates,galaxy_flux,galaxy_invvar,galaxy_z,coeff,ematrix,zvals, $
    rmatrix=rmatrix
  k_model_fluxes,coeff,galaxy_z,model_flux,ematrix=ematrix,zvals=zvals, $
    rmatrix=rmatrix
  for k=0l, nk-1 do $
    print,string(k)+' : '+string(djsig(model_flux[k,*]/galaxy_flux[k,*], $
                                       sigrej=5))
  chi2=total((model_flux-galaxy_flux)^2*galaxy_invvar,/double)
  help,chi2
  if(keyword_set(plotfluxes)) then begin
      !p.multi=[0,1,nk]
      for k=0l, nk-1 do begin
          result=moment(model_flux[k,*]/galaxy_flux[k,*])
          plot,galaxy_z,model_flux[k,*]/galaxy_flux[k,*],yra=[1.-3.*sqrt(result[1]), 1.+3.*sqrt(result[1])],psym=3,xst=1,yst=1,xra=[0.,0.5]
      endfor
  endif
endfor

; Output results

end
;------------------------------------------------------------------------------
