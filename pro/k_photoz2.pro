;+
; NAME:
;   k_photoz2.pro
;
; PURPOSE:
;   Computes the best fit template and chi2 
;   for an array of galaxy values at a given 
;   redshift
; CALLING SEQUENCE:
;
; INPUTS:
;
; OPTIONAL INPUTS:
;
; OUTPUTS:
;
; OPTIONAL INPUT/OUTPUTS:
;
; COMMENTS:
;
; EXAMPLES:
;
; BUGS:
;
; PROCEDURES CALLED:
;    k_project_filters
; REVISION HISTORY:
;   25-April-2003 Nikhil Padmanabhan (Pton), Mike Blanton (NYU)
;-
;------------------------------------------------------------------------------
pro k_photoz2, maggies, maggie_ivar, z, template_matrix, $
	       fit, chi2, $
	       filterlist = filterlist, filterpath=filterpath

;Some initial bookkeeping steps

if(NOT keyword_set(filterpath)) then $
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'
if(NOT keyword_set(filterlist)) then $
  filterlist=['sdss_u0','sdss_g0','sdss_r0','sdss_i0','sdss_z0']
foo = size(template_matrix)
n_template = (n_elements(template_matrix)/foo[1]) -1

foo = size(maggies)
n_gal = n_elements(maggies)/foo[1]

; Define some basic arrays
fit = dblarr(n_template, ngal)
chi2 = dblarr(ngal)
template_maggies = dblarr(foo[1],n_template)

; Shift lambda to the appropriate redshift

lambdas = lambda*(1.+z)

; Project the templates onto the filter curves....

template_maggies = k_project_filters,lambda,template_matrix, $
		   filterlist=filterlist, filterpath=filterpath

; Solve for the best fit using normal equations - Press et al

; Loop over all the galaxies individually
for igal = 0L, n_gal-1L do begin

; Subtract the mean spectrum from all of these
  dmaggies = maggies[*,igal] - template_maggies[*,0]
; Construct design matrix
  design = dblarr(foo[1],n_template) 
  for j = 1L,n_template do begin
   design[*,j-1] = $
	template_maggies[*,j]*sqrt(maggie_ivar[*,igal])
  endfor
; Renormalize fluxes
  dmaggies = dmaggies*sqrt(maggie_ivar[*,igal])
  dmaggies1 = transpose(design)#dmaggies
  dm  = transpose(design)#design
  fit[*,igal] = invert(dm)#dmaggies1
  chi2[igal] =  total(design#fit[*,igal])-dmaggies)^2)
endfor

; All done
end
