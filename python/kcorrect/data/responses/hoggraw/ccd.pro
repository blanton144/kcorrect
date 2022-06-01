; ccd: Output interpolated CCD response onto vector of wavelengths.
;
pro ccd, lambda, eff

    lam= [2500,3000,3500,4000,4500,5000,5500,6000,6500, $
          7000,7500,8000,8500,9000,9500,10000,10500,11000]
    eta= [0.16,0.23,0.31,0.49,0.58,0.62,0.64,0.65,0.64, $
          0.63,0.58,0.52,0.43,0.33,0.20,0.10,0.02,0.00]

    eff= ((lambda GE 2500.0) AND (lambda LE 11000.0))* $
        interpol(eta,lam,lambda)

end
;
; -----------------------------------------------------------------------------
; plotccd: Check ccd.
;
pro plotccd

    lambda= 2000.0+findgen(6000)*2.0
    ccd, lambda, eff
    plot, lambda, eff

end
