; hayes: Read and interpolate hayes data onto a vector of lambdas.
; Output in lambda f_lambda.  This extrapolates slightly, but is
; not appropriate to use in the near infrared.
;
pro hayes, lambda, lamflam

	; read in hayes data
	template= {version: 1.0, datastart: long(0), delimiter: byte(32), $
		missingvalue: 0.0, commentsymbol: "#", fieldcount: 2, fieldtypes: [4,4], $
		fieldnames: ["wavelength","relmag"], fieldlocations: [0,0], $
		fieldgroups: [0,1]}
	data= read_ascii("../hayes/hayes.txt", template=template)
    hlam= data.wavelength

	; put into flux units
    hlamflam= hlam*4.65*10.0^(-0.4*data.relmag)

    lamflam= ((lambda GE 3300.0) AND (lambda LE 10500.0))* $
        interpol(hlamflam,hlam,lambda)
    lamflam= lamflam+(lambda LT 3300.0)*hlamflam(0)
    lamflam= lamflam+(lambda GT 10500.0)*hlamflam(288)*(lambda/10500.0)^(-2.)

end
;
; -----------------------------------------------------------------------------
; plothayes: Check hayes
;
pro plothayes

    lambda= 2000.0+findgen(6000)*2.0
    hayes, lambda, lamflam
    plot, lambda, lamflam

end
