; ccdcal: calibrate a bandpass using hayes.pro, oke.pro, ccd.pro
; kwv is the water-vapor constant measured relative to typical Palomar
; conditions.  airmass is the airmass measured relative to 1.0 at
; Palomar.
;
pro ccdcal, filename, kwv, airmass

  ; read in filter information
  template= {version: 1.0, datastart: long(0), delimiter: byte(32), $
    missingvalue: 0.0, commentsymbol: "#", fieldcount: 2, fieldtypes: [4,4], $
    fieldnames: ["wavelength","transmission"], fieldlocations: [0,0], $
    fieldgroups: [0,1]}
  data= read_ascii(filename, template=template)
  lambda= data.wavelength
  dlambda= lambda-shift(lambda,1)
  dlambda(0)= 0.0
  filter= data.transmission

  ; get hayes, oke and ccd information
  hayes, lambda, hlamflam
  plot, lambda, hlamflam
  oke, lambda, kwv, airmass, okeex
  ccd, lambda, ccdeff

  ; make true throughput
  product= filter*ccdeff*10^(-0.4*okeex)
  oplot, lambda, product*max(hlamflam)/max(product)

  ; find effective wavelength
  lglameff= total(alog10(lambda)*product*dlambda)/total(product*dlambda)
  lameff= 10.0^lglameff
  oplot, [lameff,lameff],[0.0,1.0e10]
  help, lameff

  ; calibrate with nu f_nu propto nu^index spectrum
  index= 0.0
  calspec= (lambda/lameff)^(-index)
  cal= total(hlamflam*product*dlambda)/total(calspec*product*dlambda)
  calspec= calspec*cal
  oplot, lambda, calspec

  ; output in lg nu f_nu, W/m^2
  lgcal= alog10(cal)-12.0
  help, lgcal

  ; output in Jy
  Jy= cal*lameff/2.998e4
  help, Jy

  ; output in the mag of a 1 muJy source
  magmuJy= -2.5*alog10(1.0e-6/Jy)
  help, magmuJy

end
