; J. Moustakas, 2010-Jan-15, UCSD - code to convert the MMT/MEGACAM
;   filter curves to K-correct format; see the README in the parent
;   directory 
pro parse_mmt_megacam_filters

    bands = ['u','g','r','i','z']
    infilterlist = getenv('KCORRECT_DIR')+'/data/filters/mmt_megacam/'+$
      'filterthroughput.megacam'+bands+'.txt'
    outfilterlist = getenv('KCORRECT_DIR')+'/data/filters/'+$
      'mmt_megacam_'+bands+'.par'
    minwave = [3000,3500,5000,6500,7500]
    maxwave = [4000,6000,8000,9500,1.2E4]
    
    for ii = 0, n_elements(infilterlist)-1 do begin
       readcol, infilterlist[ii], skip=2, wave, resp, $
         format='F,x,x,x,x,x,x,x,x,F'
       good = where((wave gt minwave[ii]) and (wave lt maxwave[ii]),ngood)
       out = replicate({KFILTER, lambda: 0.0, pass: 0.0},ngood)
       out.lambda = wave[good]
       out.pass = resp[good]

       hdr = [$
         '# Units:',$
         '#  "lambda" is in Angstroms',$
         '#  "pass" is the contribution to the detector signal per photon',$
         '#         entering the atmosphere of Earth',$
         '#',$
         '# Bandpass Name(s): '+bands[ii],$
         '#',$
         '# Instrument: MMT/Megacam ',$
         '#',$
         '# Determined by: ',$
         '#',$
         '# Date of determination: ???',$
         '#',$
         '# Meaning of/Reason for default column: Only column available',$
         '#',$
         '# Notes:',$
         '#',$
         '# See the documentation in the MMT_MEGACAM directory',$
         '#',$
         '# These curves correspond to the total transmission, a product',$
         '# of the T(transmission)*QE(CCD Quantum efficience)*ref(Telescope',$
         '# throughput) plus atmosphere at 1.2 airmass.',$
         '#',$
         '# J. Moustakas, UCSD 2010-01-13']
       yanny_write, outfilterlist[ii], ptr_new(out), hdr=hdr
    endfor

return
end
