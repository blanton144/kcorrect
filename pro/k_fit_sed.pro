;+
; NAME:
;   k_fit_sed
;
; PURPOSE:
;   Fit SED using k_fit_coeffs and k_fit_templates iteratively
;
; CALLING SEQUENCE:
;   k_fit_sed, galaxy_maggies, galaxy_invvar, galaxy_z, templatelist, $
;     filterlist, coeff, ematrix, bmatrix, bflux, lambda, $
;     [smoothtemplate=, sublmin=, sublmax=, vmatrix=, $
;     preset_ematrix=, maxiter=, nt=, reconstruct_maggies==, plotmaggies=]
;
; INPUTS:
;   galaxy_maggies   - maggies in each band for each galaxy [N_band, N_gal]
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
;   reconstruct_maggies       - reconstructed maggies
;   plot_maggies      - set to plot pluxes as fit is done for debugging
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
pro k_fit_sed, galaxy_maggies, galaxy_invvar, galaxy_z, templatelist, filterlist, coeff, ematrix, bmatrix, bflux, lambda, smoothtemplate=smoothtemplate, sublmin=sublmin, sublmax=sublmax, vmatrix=vmatrix, preset_ematrix=preset_ematrix, maxiter=maxiter, nt=nt, reconstruct_maggies=reconstruct_maggies, plotmaggies=plotmaggies, cutlmin=cutlmin, cutlmax=cutlmax, subsmoothtemplate=subsmoothtemplate, subsmoothlimits=subsmoothlimits, useconstraint=useconstraint

; Need at least 10 parameters
if (N_params() LT 10) then begin
    klog, 'Syntax - k_fit_sed, galaxy_maggies, galaxy_invvar, galaxy_z, templatelist, $'
    klog, '    filterlist, coeff, ematrix, bmatrix, bflux, lambda, [smoothtemplate=, $'
    klog, '    sublmin=, sublmax=, vmatrix=, preset_ematrix=, maxiter=, nt=, $'
    klog, '    reconstruct_maggies=, /plotmaggies]'
    return
endif

; Set defaults
if (NOT keyword_set(maxiter)) then maxiter=100l
if (NOT keyword_set(sublmin)) then sublmin=3500.d
if (NOT keyword_set(sublmax)) then sublmax=7500.d
if (NOT keyword_set(cutlmin)) then cutlmin=2000.d
if (NOT keyword_set(cutlmax)) then cutlmax=10500.d

; Find vmatrix and bmatrix
k_load_templates,getenv('KCORRECT_DIR')+'/data/seds/'+templatelist,vmatrix, $
  lambda,cutlmin=cutlmin,cutlmax=cutlmax
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
nk=long(n_elements(galaxy_maggies))/ngalaxy
if (NOT keyword_set(nt)) then nt=nk-1l
if (NOT keyword_set(preset_ematrix)) then begin
  ematrix=dblarr(nb,nt)
  indx=lindgen(nt)
	ematrix[indx,indx]=1.d
endif
klog,'start iteration'
k_ortho_etemplates,ematrix,bflux
k_fit_coeffs,galaxy_maggies,galaxy_invvar,galaxy_z,coeff,ematrix=ematrix, $
   zvals=zvals, rmatrix=rmatrix
k_pca_etemplates,coeff,ematrix,bflux
for i=0l, maxiter-1l do begin
  k_fit_templates,galaxy_maggies,galaxy_invvar,galaxy_z,coeff,ematrix,zvals, $
    rmatrix=rmatrix
  k_reconstruct_maggies,coeff,galaxy_z,reconstruct_maggies,ematrix=ematrix, $
    zvals=zvals, rmatrix=rmatrix
  for k=0l, nk-1 do $
    print,string(k)+' : '+ $
	     string(djsig(reconstruct_maggies[k,*]/galaxy_maggies[k,*],sigrej=5))
  chi2=total((reconstruct_maggies-galaxy_maggies)^2*galaxy_invvar,/double)
  help,chi2
  if(keyword_set(plotmaggies)) then begin
      !p.multi=[0,1,nk]
      for k=0l, nk-1 do begin
          result=moment(reconstruct_maggies[k,*]/galaxy_maggies[k,*])
          plot,galaxy_z,reconstruct_maggies[k,*]/galaxy_maggies[k,*],yra=[1.-3.*sqrt(result[1]), 1.+3.*sqrt(result[1])],psym=3,xst=1,yst=1,xra=[0.,0.5]
      endfor
  endif
  k_ortho_etemplates,ematrix,bflux
	if(keyword_set(useconstraint)) then begin
		constraints_amp=1.d
	  ngals=n_elements(galaxy_z)
		scaled=coeff[1:nt-1L,*]/(replicate(1.,nt-1)#coeff[0,*])
		constraints_mean=total(scaled,2,/double)/double(ngals)
		constraints_var= 0d
		for j=0L, ngals-1L do begin 
      delta=scaled[*,j]-constraints_mean 
      constraints_var= constraints_var+delta#delta 
    endfor
		constraints_var=8.*constraints_var/double(ngals)
  endif
  k_fit_coeffs,galaxy_maggies,galaxy_invvar,galaxy_z,coeff,ematrix=ematrix, $
     zvals=zvals, rmatrix=rmatrix, constraints_amp=constraints_amp, $
     constraints_mean=constraints_mean, constraints_var=constraints_var
  k_pca_etemplates,coeff,ematrix,bflux
endfor

; Output results

end
;------------------------------------------------------------------------------
