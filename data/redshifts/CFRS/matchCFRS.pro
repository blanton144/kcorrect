;+
; NAME:
;   matchCFRS
;
; PURPOSE:
;   Matches the CFRS catalogue to the Sloan data
;
; CALLING SEQUENCE:
;   matchCFRS, rerun,outfile, size, delta=delta
;
; INPUTS:
;   rerun -- rerun to use
;   outfile -- output file
;   size -- number of input galaxies
;   delta -- matchlength in arcseconds
;
; OUTPUTS:
;   [outputfile -- FITS]
;
; EXTERNAL FUNCTIONS CALLED
;   sdss_findimage
;   sdss_readobj
;   spherematch
;
; COMMENTS:
;   1. Simply uses only the best revisions that sdss_findimage returns
;   2. Requires that you specify the number of CFRS galaxies (still
;	needs to be fixed)
; REVISION HISTORY:
;   N. Padmanabhan 2003-May-01
;
;----------------------------------------------------------------------

pro matchCFRS, rerun, outfile, size, delta=delta

; Read in the CFRS data
  data = dblarr(3,size)
  openr, lun, "CFRS", /get_lun
  readf, lun, data
  free_lun, lun

; Trim out stars and unidentified sources
; these have z=0.0,9.9 ( we use >5.0)
; Also, define the structure

  index = where((data[2,*] GT 0.00) AND (data[2,*] LT 5.0))
  nel  = n_elements(index)
  res = {ra:0.0d0, dec:0.0d0, z:0.0d0, $
    matched:0, run:0L, petroflux:dblarr(5), petroflux_ivar:dblarr(5), $
    modelflux:dblarr(5), modelflux_ivar:dblarr(5), extinction:dblarr(5)}
  struct = replicate(res,nel)
; Strange IDLism hack
  tmp = data[0,index]
  struct[*].ra = tmp[*]
  tmp = data[1,index]
  struct[*].dec = tmp[*]
  tmp = data[2,index]
  struct[*].z = tmp[*]
  data = ""

; Get the runs that fit
  res = ""
  res = sdss_findimage(struct[*].ra, struct[*].dec, /best, rerun=rerun)

; Now actually read the data
  found = where(res[*].run GT 0)

; Must only extract the unique elements
; This appears to be a non-trivial statement
  runndx = uniq(res[found].run, sort(res[found].run))
  count = n_elements(runndx)
  for i = 0L,count-1 do begin
    run = (res[found].run)[runndx[i]]
    tfield = (res[found].field)[where(res[found].run eq run)]  
    field = tfield[uniq(tfield,sort(tfield))]
    for camcol = 1,6 do begin
     sdss = sdss_readobj(run,camcol,field, $
        rerun = rerun, select_tags=['RA','DEC','PETROFLUX*','MODELFLUX*','EXTINCTION'])
  ; Did this return anything?
     if (keyword_set(sdss) EQ 0) then CONTINUE
  ; Now do spherematches
     spherematch, struct[*].ra, struct[*].dec, sdss[*].ra, sdss[*].dec, $
       delta/3600.0d0, match1, match2, distance12
     nmatch = n_elements(match1)
     if (nmatch EQ 0) then CONTINUE 
     struct[match1].matched = 1   ; Marked as matched
     struct[match1].petroflux = sdss[match2].petroflux
     struct[match1].petroflux_ivar = sdss[match2].petroflux_ivar
     struct[match1].modelflux = sdss[match2].modelflux
     struct[match1].modelflux_ivar = sdss[match2].modelflux_ivar
     struct[match1].extinction = sdss[match2].extinction
    endfor
   endfor

; Copy data

  print, "Writing out file : ",outfile
  mwrfits, struct, outfile, /Create

  index = where(struct[*].matched EQ 1)
  print, n_elements(index)
  
end









