;+
;NAME:
;  k_read_dale.pro
;PURPOSE:
;  reads the Dale SEDs
;CALLING SEQUENCE:
;  dale=k_read_dale()
;OUTPUTS:
;  dale - structure with 
;             .LAMBDA - wavelength in angstroms
;             .FLUX - flux 
;BUGS:
;  What units are fluxes in?
;  Untested?
;REVISION HISTORY:
;  2004-Dec-20  started by Hogg
;_
;------------------------------------------------------------
function k_read_dale
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
dale.lambda= transpose(blah.field01[0,*])*1D3
dale.flux= blah.field01[1:nsed,*]
return, dale
end
