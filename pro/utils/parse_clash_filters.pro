; J. Moustakas, 2011-Apr-07, UCSD - code to convert the HST ACS/WFC
; and WFC3 filter response functions used by the CLASH collaboration
; to K-correct format; see the README in the data/filters directory
; for more details

pro parse_clash_filters

    filtdir = getenv('KCORRECT_DIR')+'/data/filters/'
    infilterlist = file_search(filtdir+'clash/*.res',count=nfilt)

    root = repstr(repstr(repstr(repstr(repstr($
      file_basename(infilterlist),'HST_'),'.res',''),'WFC_'),'IR_'),'UVIS_')
    instr = repstr(repstr(file_basename(infilterlist),'HST_'),'.res','')
    bands = repstr(strmid(repstr(repstr(file_basename(infilterlist),$
      'HST_'),'.res',''),5,6,/rever),'_','')
    
    outfilterlist = filtdir+'clash_'+strlowcase(root)+'.par'
    
    for ii = 0, n_elements(infilterlist)-1 do begin
       readcol, infilterlist[ii], wave, resp, format='F,F', /silent

       keep = where(resp gt 1D-4,nwave)
;      nwave = n_elements(wave)
       out = replicate({KFILTER, lambda: 0.0, pass: 0.0},nwave)
       out.lambda = wave[keep]
       out.pass = resp[keep]
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
         '# Instrument: '+instr[ii],$
         '#',$
         '# Determined by: ',$
         '#',$
         '# Date of determination: 2010-May-02',$
         '#',$
         '# Meaning of/Reason for default column: Only column available',$
         '#',$
         '# Notes:',$
         '#',$
         '#',$
         '# Obtained from Leonidas Moustakas on 2011-Apr-02 for the CLASH team.',$
         '#',$
         '# J. Moustakas, UCSD 2011-Apr-07']
       yanny_write, outfilterlist[ii], ptr_new(out), hdr=hdr
    endfor

return
end
