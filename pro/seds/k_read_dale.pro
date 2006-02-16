;+
;NAME:
;  k_read_dale.pro
;PURPOSE:
;  reads the Dale SEDs
;CALLING SEQUENCE:
;  dale=k_read_dale( [alpha= ] )
;OUTPUTS:
;  dale - structure with 
;             .LAMBDA - wavelength in angstroms
;             .FLUX[64] - flux for each value of ALPHA
;  alpha - values of alpha
;BUGS:
;  What units are fluxes in?
;  Untested?
;REVISION HISTORY:
;  2004-Dec-20  started by Hogg
;_
;------------------------------------------------------------
function k_read_dale, alpha=alpha

dale_file= getenv('KCORRECT_DIR')+'/data/seds/dale/spectra.dat'
blah= read_ascii(dale_file)
tmp= size(blah.field01,/dimens)
nsed= tmp[0]-1
nlam= tmp[1]
dale= {dale, $
       lambda     : 0D,          $
       flux       : dblarr(nsed) $
      }
dale= replicate(dale,nlam)
dale.lambda= transpose(blah.field01[0,*])*1D4
for i=0, nlam-1 do $
  dale[i].flux= 10.^(blah.field01[1:nsed,i])/dale[i].lambda

alphastr=read_ascii(getenv('KCORRECT_DIR')+'/data/seds/dale/alpha.dat')
alpha=transpose(alphastr.field1[0,*])

return, dale
end
