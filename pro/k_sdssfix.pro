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
pro k_sdssfix,mags,magserr,galaxy_z=galaxy_z, errlimit=errlimit, maglimit=maglimit, errband=errband, defaultcoeff=defaultcoeff, defaultmags=defaultmags, defaultz=defaultz, version=version, vpath=vpath, rmatrix=rmatrix, zvals=zvals, ematrix=ematrix, filterpath=filterpath, defaultzlimits=defaultzlimits

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
if(n_elements(errlimit) eq 0) then errlimit=[0.8,0.8,0.8,0.8,0.8]
if(n_elements(errband) eq 0) then errband=[0.05,0.02,0.02,0.02,0.03]
if(n_elements(maglimit) eq 0) then maglimit=[23.0,23.0,23.0,23.0,23.0]

; Make sure we are working in magnitudes; negative or zero maggies
; are just errors so we set them to the limits
if(keyword_set(maggies)) then begin 
  for k=0, nk-1l do begin
    posindx=where(mags[k,*] gt 0.,poscount)
    negindx=where(mags[k,*] le 0.,negcount)
    if(negcount gt 0) then begin
	     magserr[k,negindx]=errlimit[k]
	     mags[k,negindx]=maglimit[k]
	  endif
    if(poscount gt 0) then begin
		   if(keyword_set(invvar)) then $
				     magserr[k,posindx]=sqrt(1./magserr[k,posindx])
			 magserr[k,posindx]=2.5d*alog(0.1d)*magserr[k,posindx]/mags[k,posindx]
			 mags[k,posindx]=-2.5*alog10(mags[k,posindx])
		endif
  endfor
endif

; First, fix the errors 
for k=0, nk-1 do magserr[k,*]=sqrt(magserr[k,*]^2+errband[k]^2)

; Second, determine which galaxies and bands need fixing; absolute magnitudes
; are there to take care of PHOTO's -9999 stuff
badgalaxy=lonarr(ngalaxy)
badband=lonarr(nk,ngalaxy)
for k=0, nk-1 do begin
  badindx=where(abs(mags[k,*]) ge maglimit[k] or $
								abs(magserr[k,*]) ge errlimit[k],count)
	if(count gt 0) then begin
	   badgalaxy[badindx]=1l
	   badband[k,badindx]=1l
  endif
endfor
badindx=where(badgalaxy gt 0, count)
if(count eq 0) then return

; Third, set all the errors in all of the bad bands to the appropriate 
; limiting error
for k=0, nk-1l do begin
  badindx=where(badband[k,*] gt 0,count)
	if(count gt 0) then begin
     magserr[k,badindx]=errlimit[k]
  endif
endfor 

; Fourth, set the "reasonable" coeffs -- the coefficients for a typical 
; galaxy SED (based on its z=0.1 SDSS colors); we wait to do this until
; now just in case there are no bad objects
if(n_elements(defaultz) eq 0) then defaultz=0.1
if(n_elements(defaultmags) eq 0) then $
  defaultmags=[-10.155465,-11.563362,-12.356318,-12.749081,-13.017499]
if(n_elements(defaultcoeff) eq 0) then begin
  defaultmaggies=10.^(-0.4*(defaultmags))
  defaultinvvar=1./(defaultmaggies*0.02)^2
	if(n_elements(rmatrix) gt 0 AND n_elements(zvals) gt 0 AND $
		 n_elements(ematrix) gt 0) then begin
	   k_fit_coeffs,defaultmaggies,defaultinvvar,defaultz,defaultcoeff, $
	      filterpath=filterpath,rmatrix=rmatrix,zvals=zvals, $
	      ematrix=ematrix
	endif else begin
	   k_fit_coeffs,defaultmaggies,defaultinvvar,defaultz,defaultcoeff, $
	    version=version,vpath=vpath,filterpath=filterpath,rmatrix=rmatrix, $
      zvals=zvals,ematrix=ematrix
	endelse
endif
usedefaultcoeff=dblarr(n_elements(defaultcoeff),n_elements(to_z))
for t=0, n_elements(defaultcoeff)-1 do begin
  usedefaultcoeff[t,*]=defaultcoeff[t]
endfor 
if(n_elements(rmatrix) gt 0 AND n_elements(zvals) gt 0 AND $
	 n_elements(ematrix) gt 0) then begin
  k_reconstruct_maggies,usedefaultcoeff,to_z,reconstruct_maggies, $
    filterpath=filterpath,rmatrix=rmatrix,zvals=zvals, $
    ematrix=ematrix
endif else begin
  k_reconstruct_maggies,usedefaultcoeff,to_z,reconstruct_maggies, $
    version=version,vpath=vpath,filterpath=filterpath,rmatrix=rmatrix, $
    zvals=zvals,ematrix=ematrix
endelse
fixmags=dblarr(nk,n_elements(to_z))
for k=0, nk-1 do begin
  negindx=where(reconstruct_maggies[k,*] le 0.d,negcount)
	posindx=where(reconstruct_maggies[k,*] gt 0.d,poscount)
	if(negcount gt 0) then $
    fixmags[k,negindx]=maglimit[k]
	if(poscount gt 0) then $
    fixmags[k,posindx]=-2.5*alog10(reconstruct_maggies[k,*])
endfor  
if(n_elements(to_z) eq 1) then begin
  tmpfixmags=fixmags
  fixmags=dblarr(nk,ngalaxy)
	for k=0, nk-1 do fixmags[k,*]=tmpfixmags[k]
endif

; Fifth, fix the magnitudes, as follows: 
;   for each band from u to i, determine the nearest, good, redder band for 
;     each bad object. if such a band exists, fix based on it. Call that band
;     good now
;   for each band from g to z, determine the nearest, good, bluer band for 
;     each bad object
;   the only remaining objects should be all bad bands. Give these the
;     defaultmags, all relative to the r-band
for k=0, nk-2 do begin
  for kp=k+1, nk-1 do begin
    fixindx=where(badband[k,*] gt 0 and badband[kp,*] eq 0, count)
	  if(count gt 0) then begin
		  mags[k,fixindx]=fixmags[k,fixindx]-fixmags[kp,fixindx]+mags[kp,fixindx]
      badband[k,fixindx]=0l
		endif
	endfor
endfor
for k=nk-1, 1, -1 do begin
  for kp=k-1, 0, -1 do begin
    fixindx=where(badband[k,*] gt 0 and badband[kp,*] eq 0, count)
	  if(count gt 0) then begin
		  mags[k,fixindx]=fixmags[k,fixindx]-fixmags[kp,fixindx]+mags[kp,fixindx]
      badband[k,fixindx]=0l
		endif
	endfor
endfor
stillbadfixr=where(total(badband,1) gt 0 and abs(mags[2,*]) gt maglimit[2], $
                   count)
if(count gt 0) then $
  mags[2,stillbadfixr]=maglimit[2]
stillbad=where(total(badband,1) gt 0,count)
if(count gt 0) then begin
  for k=0, nk-1 do begin
    if(k ne 2) then $
      mags[k,stillbad]=fixmags[k,stillbad]-fixmags[2,stillbad]+mags[2,stillbad]
  endfor
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
