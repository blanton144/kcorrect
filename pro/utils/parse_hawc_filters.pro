; J. Moustakas, 2015-May-19 - code to convert the Sofia/HAWC filter response
; functions to K-correct format; see the README in the data/filters directory
; for more details

pro parse_hawc_filters

    filtdir = getenv('KCORRECT_DIR')+'/data/filters/'
    infilterlist = file_search(filtdir+'hawc/band?.txt',count=nfilt)
    
    bands = 'Band_'+strtrim(indgen(4)+1,2)
    outfilterlist = filtdir+'hawc_band_'+strtrim(indgen(4)+1,2)+'.par'
    
    for ii = 0, n_elements(infilterlist)-1 do begin
       readcol, infilterlist[ii], wave, resp, format='F,X,X,F', /silent, skip=2

       nwave = n_elements(wave)
       out = replicate({KFILTER, lambda: 0.0, pass: 0.0},nwave)
       out.lambda = wave*1D4
       out.pass = resp

;      plot, out.lambda, out.pass, ysty=3, xsty=3
;      cc = get_kbrd(1)
       
       hdr = [$
         '# Units:',$
         '#  "lambda" is in Angstroms',$
         '#  "pass" is the contribution to the detector signal per photon',$
         '#         entering the atmosphere of Earth',$
         '#',$
         '# Bandpass Name(s): '+bands[ii],$
         '#',$
         '# Instrument: Sofia/HAWC',$
         '#',$
         '# Determined by: HAWC Instrument Team',$
         '#',$
         '# Date of determination: Unknown',$
         '#',$
         '# Meaning of/Reason for default column: Only column available',$
         '#',$
         '# Notes:',$
         '#',$
         '#',$
         '# Obtained from Ravi Sankrit at the Sofia helpdesk on 2015-May-19 on ',$
         '# behalf of the HAWC Instrument Team.  Transmission includes the effect',$
         '# of the atmosphere but it is not clear if it includes the CCD QE.',$
         '#',$
         '# J. Moustakas, 2015-May-19, Siena College']
       yanny_write, outfilterlist[ii], ptr_new(out), hdr=hdr
    endfor

return
end
