;+
; NAME:
;   k_viewer_image_fit
; PURPOSE:
;   Run viewer image fit in a directory
; CALLING SEQUENCE:
; INPUTS:
; OPTIONAL INPUTS:
; KEYWORDS:
; OUTPUTS:
; OPTIONAL INPUT/OUTPUTS:
; COMMENTS:
; EXAMPLES:
; BUGS:
; PROCEDURES CALLED:
; REVISION HISTORY:
;   25-jan-2002  WRitten by Mike Blanton, NYU
;-
;------------------------------------------------------------------------------
pro k_run_vif

filenames=findfile('image*.fit')

for i=0L, n_elements(filenames)-1L do begin
    k_viewer_image_fit,filenames[i],red,blue
    name=(stregex(filenames[i],'image(.*)\.fit',/subexpr,/extract))[1]
    redbluefile='redblue'+name+'.fits'
    mwrfits,red,redbluefile,/create
    mwrfits,blue,redbluefile
    nx=(size(red,/dimensions))[0]
    ny=(size(red,/dimensions))[1]
    lohi=dblarr(2)
    lohi[1]=max([red[nx/2-10:nx/2+10,ny/2-10:ny/2+10], $
                 blue[nx/2-10:nx/2+10,ny/2-10:ny/2+10]])*0.85
    lohi[0]=-0.1*lohi[1]
    if(i eq 0) then $
      greyscale,red,0.396,'arcsec',lohi[1],lohi[0],'run_vif_qa.ps', $
      title=filenames[i]+' (red)',/startps 
    hogg_greyscale_plot,red,pixscale=0.396,scalename='arcsec',lo=lohi[1], $
      hi=lohi[0],title=filenames[i]+' (red)' 
    hogg_greyscale_plot,blue,pixscale=0.396,scalename='arcsec',lo=lohi[1], $
      hi=lohi[0],title=filenames[i]+' (blue)'
endfor
end_print

end
;------------------------------------------------------------------------------
