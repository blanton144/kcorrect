pro k_create_plots

;run_fit_coeffs,'default',savfile='all.default.sav',/evenz, $
	;modelzlim=0.25
;run_fit_coeffs,'default',savfile='main.default.sav',primtargetmask=64, $
;	modelzlim=0.5
;run_fit_coeffs,'default',savfile='lrg.default.sav',primtargetmask=32, $
	;modelzlim=0.01

k_speck,'main.default.sav','all.default.specK.0.1.sav',subsample=500, to_z=0.1
$touch gotpast
k_speck,'lrg.default.sav','all.default.specK.0.3.sav',subsample=500, to_z=0.3

;k_coeffdist_plot,'all.default.sav',psfile='k_coeffdist_plot.ps', $
;  subsample=3,basecoeff=1
k_colors_plot,'main.default.sav',psfile='main_colors_plot.ps', $
  zrange=[0.,0.23], to_z=0.1,lumlim=[-21.5,-21.2], $
  colorlimits=[[0.5,2.7],[0.23,1.13],[0.25,0.5],[0.01,0.52]], $
  subsample=3, zsplit=[0.05,0.1,0.17]
k_colors_plot,'lrg.default.sav',psfile='lrg_colors_plot.ps', $
  zrange=[0.1,0.49], to_z=0.3,lumlim=[-22.8,-22.5], $
  nprimtargetmask=67108864, $
  colorlimits=[[0.7,2.9],[0.84,1.74],[0.43,0.68],[0.09,0.59]], $
  zsplit=[0.2,0.3,0.4]
k_speck_plot,'all.default.specK.0.1.sav',to_z=0.1,nsig=3.4, $
  psfile='k_speck_plot.fiber.0.1.ps',usefiber='main.default.sav', $
  zrange=[0.02,0.21]
k_speck_plot,'all.default.specK.0.3.sav',to_z=0.3,nsig=1.5, $
  psfile='k_speck_plot.fiber.0.3.ps',usefiber='lrg.default.sav', $
  zrange=[0.1,0.5],lumlim=[-24.0,-21.5]
k_kcorrect_plot,'all.default.sav',psfile='k_kcorrect_plot.ps',subsample=4, $
  to_z=0.3
k_kcorrect_grdiff_plot,'all.default.sav',psfile='k_kcorrect_grdiff_plot.ps', $
  subsample=4,to_z=0.3
k_speck_plot,'all.default.specK.0.1.sav',to_z=0.1,ylimits=[-0.16,0.16], $
  psfile='k_speck_plot.0.1.ps',zrange=[0.02,0.21],primtargetmask=64, $
  lumlim=[-23.,-19.]
k_speck_plot,'all.default.specK.0.3.sav',to_z=0.3,ylimits=[-0.29,0.29], $
  psfile='k_speck_plot.0.3.ps',zrange=[0.1,0.5]
k_speck_plot,'all.default.specK.0.1.sav',to_z=0.1,ylimits=[-0.16,0.16], $
  psfile='k_speck_plot.fitfib.0.1.ps',zrange=[0.02,0.21],primtargetmask=64, $
  lumlim=[-23.,-19.],/fitfib
k_speck_plot,'all.default.specK.0.3.sav',to_z=0.3,ylimits=[-0.29,0.29], $
  psfile='k_speck_plot.fitfib.0.3.ps',zrange=[0.1,0.5],/fitfib

;$cp -f k_kcorrect_plot.ps main_colors_plot.ps lrg_colors_plot.ps k_coeffdist_plot.ps k_speck_plot.0.1.ps k_speck_plot.0.3.ps k_speck_plot.fitfib.0.1.ps k_speck_plot.fitfib.0.3.ps $KCORRECT_DIR/docs/paper

end
