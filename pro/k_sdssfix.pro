;+
; NAME:
;   k_sdssfix
;
; PURPOSE:
;   Take a set of SDSS magnitudes or maggies, and errors, and "fix"
;   them, in the sense that for "bad" measurements or errors you 
;   assign values which are not absurd, with large error bars. 
;  
;   1. Fixes the errors to be reasonable
;   2. Fixes "bad" colors
;
; CALLING SEQUENCE:
;   k_sdssfix, mags, magserr [, z=z, errlimit=, maglimit=, errband=, $
;      defaultcoeff=, version=, vpath= ]
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
;
; REVISION HISTORY:
;   07-Feb-2002  Written by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_sdssfix,mags,magserr,galaxy_z=galaxy_z, errlimit=errlimit, maglimit=maglimit, errband=errband, defaultcoeff=defaultcoeff, defaultmags=defaultmags, defaultz=defaultz, version=version, vpath=vpath, rmatrix=rmatrix, zvals=zvals, ematrix=ematrix, filterpath=filterpath, defaultzlimits=defaultzlimits, maggies=maggies, invvar=invvar, errorsonly=errorsonly

nk=5
ngalaxy=n_elements(mags)/nk

if(nk ne 5) then begin
  klog, 'You are supposed to be using 5-band SDSS observations!'
  return 
endif

; Set default parameters
if(NOT keyword_set(filterpath)) then $
  filterpath=getenv('KCORRECT_DIR')+'/data/filters'
if(NOT keyword_set(vpath)) then $
  vpath=getenv('KCORRECT_DIR')+'/data/etemplates'
if(NOT keyword_set(version)) then $
  version='default'
if(n_elements(galaxy_z) eq 0) then begin 
   to_z=[0.1]
endif else begin
   to_z=[galaxy_z]
endelse
if(n_elements(defaultzlimits) lt 2) then defaultzlimits=[0.,0.7]
outzlimits=where(to_z lt defaultzlimits[0],count)
if(count gt 0) then to_z[outzlimits]=defaultzlimits[0]
outzlimits=where(to_z gt defaultzlimits[1],count)
if(count gt 0) then to_z[outzlimits]=defaultzlimits[1]
if(n_elements(errlimit) eq 0) then errlimit=3.
if(n_elements(largeerrlimit) eq 0) then largeerrlimit=20.0
if(n_elements(errband) eq 0) then errband=[0.05,0.02,0.02,0.02,0.03]
if(n_elements(largemaglimit) eq 0) then largemaglimit=50.0
if(n_elements(maglimit) eq 0) then maglimit=26.

; Make sure we are working in magnitudes; negative or zero maggies
; are just errors so we set them to the limits
if(keyword_set(maggies)) then begin 
  for k=0, nk-1l do begin
    posindx=where(mags[k,*] gt 0.,poscount)
    negindx=where(mags[k,*] le 0.,negcount)
    if(negcount gt 0) then begin
	     magserr[k,negindx]=errlimit
	     mags[k,negindx]=maglimit
	  endif
    if(poscount gt 0) then begin
		   if(keyword_set(invvar)) then $
				     magserr[k,posindx]=sqrt(1./magserr[k,posindx])
			 magserr[k,posindx]=2.5d/alog(10.d)*magserr[k,posindx]/mags[k,posindx]
			 mags[k,posindx]=-2.5*alog10(mags[k,posindx])
		endif
  endfor
endif

; First, fix the errors 
for k=0, nk-1 do magserr[k,*]=sqrt(magserr[k,*]^2+errband[k]^2)

; Second, determine which galaxies and bands need fixing because they have 
; wrong errors; absolute magnitudes are there to take care of PHOTO's 
; -9999 stuff
badindx=where(abs(magserr) ge largeerrlimit,count)
if(count gt 0) then begin
  magserr[badindx]=largeerrlimit
endif

; Third, determine which need fixing because of insane magnitude values
badindx=where(abs(mags) ge largemaglimit,count)
if(count gt 0) then begin
  mags[badindx]=maglimit
  magserr[badindx]=errlimit
endif

; Fourth, determine which need fixing because of faint magnitude values
badindx=where(abs(mags) ge maglimit,count)
if(count gt 0) then begin
  magserr[badindx]=sqrt(magserr[badindx]^2+errlimit^2)
endif

; Finally, convert back to maggies if desired
if(keyword_set(maggies)) then begin 
    mags=10.d^(-0.4d*mags)
    magserr=mags*0.4d*alog(10.d)*magserr
    if(keyword_set(invvar)) then begin 
		    magserr=1.d/magserr^2
    endif
endif

end
;------------------------------------------------------------------------------
