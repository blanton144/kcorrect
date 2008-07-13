pro flamingos2filters
; J. Moustakas, 2008-Jul-11, NYU
    
    filterpath = getenv('KCORRECT_DIR')+'/data/filters/'
    datapath = filterpath+'/flamingos/'

    hdr0='# Units:'
    hdr0 = [hdr0, '#  "lambda" is in Angstroms']
    hdr0 = [hdr0, '#  "pass" is the contribution to the detector signal per photon']
    hdr0 = [hdr0, '#']
    hdr0 = [hdr0, '# Bandpass Name(s): ']
    hdr0 = [hdr0, '#']
    hdr0 = [hdr0, '# Instrument: FLAMINGOS']
    hdr0 = [hdr0, '#']
    hdr0 = [hdr0, '# Determined by: ???']
    hdr0 = [hdr0, '#']
    hdr0 = [hdr0, '# Date of determination: July 2003(?)']
    hdr0 = [hdr0, '#']
    hdr0 = [hdr0, '# Notes:']
    hdr0 = [hdr0, '#']
    hdr0 = [hdr0, '# Retrieved from http://flamingos.astro.ufl.edu/Filter_Info/index.html']
    hdr0 = [hdr0, '# Filter curve (measured within the dewar) has been convolved with the']
    hdr0 = [hdr0, '# atmosphere (AM=1.0), but does not include the detector QE or telescope']
    hdr0 = [hdr0, '# throughput.']
    hdr0 = [hdr0, '#']
    hdr0 = [hdr0, '# J. Moustakas, 2008-Jul-12, NYU']

    kfilter1 = {lambda:0.D, pass:0.D}

    bandpass = ['J','H','Ks']
    outfile = filterpath+'flamingos_'+bandpass+'.par'
    infile = datapath+'FLAMINGOS.BARR.'+['J.MAN240','H.MAN109',$
      'Ks.MAN306A']+'.ColdWitness.txt'

    for ii = 0L, n_elements(infile)-1L do begin

       readcol, infile[ii], wave, resp, format='F,F', /silent
       resp = resp/100.0 & wave = wave*10.0
       resp = resp*atmo_transmission(wave) ; convolve with the atmosphere

       hdr = hdr0
       here = where(hdr eq '# Bandpass Name(s): ')
       hdr[here] = '# Bandpass Names(s): '+bandpass[ii]
       kfilter = replicate(kfilter1,n_elements(resp))
       kfilter.lambda = wave
       kfilter.pass = resp
       yanny_write, outfile[ii], ptr_new(kfilter), $
         hdr=hdr, stnames='KFILTER'

    endfor
       
return
end
    
