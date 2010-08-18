; J. Moustakas, 2010-Aug-16, UCSD - code to convert the optical
; (Ugriz) filter curves used by the SWIRE team to image the CDFS field
; (with the CTIO/MOSAIC-II instrument) to K-correct format; see the
; README in the parent directory for more details

pro parse_cdfs_swire_filters, debug=debug

    bands = ['U','g','r','i','z']
    infilterlist = getenv('KCORRECT_DIR')+'/data/filters/cdfs_swire/'+$
      bands+'.dat'
    outfilterlist = getenv('KCORRECT_DIR')+'/data/filters/'+$
      'cdfs_swire_'+bands+'.par'

    readcol, getenv('KCORRECT_DIR')+'/data/filters/cdfs_swire/qe.dat', $
      qewave, qe, skip=1, format='F,F', /silent
    qe = qe/100.0
    
    for ii = 0, n_elements(infilterlist)-1 do begin
       readcol, infilterlist[ii], skip=1, wave, resp, $
         format='F,F'
       nwave = n_elements(wave)
       out = replicate({KFILTER, lambda: 0.0, pass: 0.0},nwave)
       out.lambda = wave
       out.pass = (resp>0)/100.0

       linterp, qewave, qe, out.lambda, qe1, missing=0.0 ; detector QE
       oke, out.lambda, 1.0, 1.0, ext ; atmosphere

       if keyword_set(debug) then begin
          djs_plot, out.lambda, out.pass, xsty=3, ysty=3
          djs_oplot, qewave, qe, line=5, psym=-6
          djs_oplot, out.lambda, out.pass*qe1, color='red'
          djs_oplot, out.lambda, out.pass*qe1*10^(-0.4*ext), color='green'
          cc = get_kbrd(1)
       endif
       out.pass = out.pass*qe1*10.0^(-0.4*ext)
       
       hdr = [$
         '# Units:',$
         '#  "lambda" is in Angstroms',$
         '#  "pass" is the contribution to the detector signal per photon',$
         '#         entering the atmosphere of Earth',$
         '#',$
         '# Bandpass Name(s): '+bands[ii],$
         '#',$
         '# Instrument: CTIO/MOSAIC-II ',$
         '#',$
         '# Determined by: ',$
         '#',$
         '# Date of determination: ???',$
         '#',$
         '# Meaning of/Reason for default column: Only column available',$
         '#',$
         '# Notes:',$
         '#',$
         '# See the documentation in the CDFS_SWIRE directory',$
         '#',$
         '# These curves correspond to the total transmission, a product',$
         '# of the T(transmission)*QE(CCD Quantum efficience)*ref(Telescope',$
         '# throughput) plus atmosphere at 1.0 airmass.',$
         '#',$
         '# J. Moustakas, UCSD 2010-08-16']
       yanny_write, outfilterlist[ii], ptr_new(out), hdr=hdr
    endfor

return
end
