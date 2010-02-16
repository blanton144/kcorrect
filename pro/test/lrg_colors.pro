pro lrg_colors

common com_sdss_kcorrect, rmatrix, zvals, band_shift, vname, ermatrix

lowz=lowz_read(sample='dr6')
indx= shuffle_indx(n_elements(lowz), num_sub=600)
lowz=lowz[indx]

kc=sdss_kcorrect(lowz.zdist, cal=lowz, coeffs=coeffs, mass=mass)

k_print, filename='lrg_colors.ps'

djs_plot, [0], [0], xra=[-0.1, 2.0], yra=[-0.1, 1.5]

redshift=0.1+0.6*(0.5+findgen(100))/100.
imag0= 19.8
hogg_usersym, 12, /fill
for i=0L, n_elements(lowz)-1L do begin
    tmp_coeffs=coeffs[*,i]#replicate(1., n_elements(redshift))
    k_reconstruct_maggies, tmp_coeffs, redshift, rmaggies, $
      rmatrix=rmatrix, zvals=zvals
    imag= -2.5*alog10(rmaggies[3,*])
    dmag= lf_distmod(redshift)-(lf_distmod(lowz[i].zdist))[0] + $
          (imag-imag0)
    tmp_mass= mass[i]*10.^(0.4*dmag)
    gmr= -2.5*alog10(rmaggies[1,*]/rmaggies[2,*])
    rmi= -2.5*alog10(rmaggies[2,*]/rmaggies[3,*])

    ii=where(tmp_mass lt 3.e+9, nii)
    if(nii gt 0) then $
      djs_oplot, gmr[ii], rmi[ii], psym=8, symsize=0.35, color='blue'

    ii=where(tmp_mass gt 3.e+9 and tmp_mass lt 1.e+10, nii)
    if(nii gt 0) then $
      djs_oplot, gmr[ii], rmi[ii], psym=8, symsize=0.35, color='cyan'
    
    ii=where(tmp_mass gt 1.e+10 and tmp_mass lt 3.e+10, nii)
    if(nii gt 0) then $
      djs_oplot, gmr[ii], rmi[ii], psym=8, symsize=0.30, color='green'

    ii=where(tmp_mass gt 3.e+10 and tmp_mass lt 6.e+10, nii)
    if(nii gt 0) then $
      djs_oplot, gmr[ii], rmi[ii], psym=8, symsize=0.25, color='yellow'

    ii=where(tmp_mass gt 6.e+10 and tmp_mass lt 1.e+11, nii)
    if(nii gt 0) then $
      djs_oplot, gmr[ii], rmi[ii], psym=8, symsize=0.20, color='orange'

    ii=where(tmp_mass gt 1.e+11, nii)
    if(nii gt 0) then $
      djs_oplot, gmr[ii], rmi[ii], psym=8, symsize=0.15, color='red'

endfor

xst= !X.CRANGE[0]+0.03*(!X.CRANGE[1]-!X.CRANGE[0])
yst= !Y.CRANGE[0]+0.95*(!Y.CRANGE[1]-!Y.CRANGE[0])
djs_xyouts, xst, yst, textoidl('M>1\times10^{11} h^{-2} M_\odot'), $
            color='red'

yst= !Y.CRANGE[0]+0.91*(!Y.CRANGE[1]-!Y.CRANGE[0])
djs_xyouts, xst, yst, textoidl('6\times10^{10} h^{-2} M_\odot<M<'+ $
                               '1\times10^{11} h^{-2} M_\odot'), $
            color='orange'

yst= !Y.CRANGE[0]+0.87*(!Y.CRANGE[1]-!Y.CRANGE[0])
djs_xyouts, xst, yst, textoidl('3\times10^{10} h^{-2} M_\odot<M<'+ $
                               '6\times10^{10} h^{-2} M_\odot'), $
            color='yellow'

yst= !Y.CRANGE[0]+0.83*(!Y.CRANGE[1]-!Y.CRANGE[0])
djs_xyouts, xst, yst, textoidl('1\times10^{10} h^{-2} M_\odot<M<'+ $
                               '3\times10^{10} h^{-2} M_\odot'), $
            color='green'

yst= !Y.CRANGE[0]+0.79*(!Y.CRANGE[1]-!Y.CRANGE[0])
djs_xyouts, xst, yst, textoidl('3\times10^{9} h^{-2} M_\odot<M<'+ $
                               '1\times10^{10} h^{-2} M_\odot'), $
            color='cyan'

yst= !Y.CRANGE[0]+0.75*(!Y.CRANGE[1]-!Y.CRANGE[0])
djs_xyouts, xst, yst, textoidl('M<3\times10^{9} h^{-2} M_\odot'), $
            color='blue'

k_end_print

end
