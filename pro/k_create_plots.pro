pro k_create_plots,version,vpath, addgrgap=addgrgap, sdssfix=sdssfix, vconstraint=vconstraint,nsp=nsp,dots=dots

if(NOT keyword_set(version)) then version='default'
if(NOT keyword_set(vconstraint)) then vconstraint=1
if(NOT keyword_set(sdssfix)) then sdssfix=1

run_fit_coeffs,version,savfile='all.'+version+'.sav',/evenz, $
	modelzlim=0.25,spfile='/globaldata/howdy3/sdss/spectro/spAll.fits', $
  outpts='all.'+version+'.pts',vpath=vpath, addgrgap=addgrgap, $
  sdssfix=sdssfix, vconstraint=vconstraint, nsp=nsp
run_fit_coeffs,version,savfile='main.'+version+'.sav',primtargetmask=64, $
	modelzlim=0.5,spfile='/globaldata/howdy3/sdss/spectro/spAll.fits', $
  outpts='main.'+version+'.pts', vpath=vpath, addgrgap=addgrgap, $
  sdssfix=sdssfix, vconstraint=vconstraint, nsp=nsp
run_fit_coeffs,version,savfile='lrg.'+version+'.sav',primtargetmask=32, $
	modelzlim=0.01,spfile='/globaldata/howdy3/sdss/spectro/spAll.fits', $
  outpts='lrg.'+version+'.pts', vpath=vpath, addgrgap=addgrgap, $
  sdssfix=sdssfix, vconstraint=vconstraint, nsp=nsp

k_colors_plot,'lrg.'+version+'.sav',psfile='lrg_colors_plot.z0.2.ps', $
  zrange=[0.1,0.49], to_z=0.2,lumlim=[-22.8,-22.5], $
  nprimtargetmask=67108864, $
  colorlimits=[[0.7,2.9],[0.84,1.74],[0.43,0.68],[0.09,0.59]], $
  zsplit=[0.2,0.3,0.4], version=version, vpath=vpath, addgrgap=addgrgap, $
  sdssfix=sdssfix, vconstraint=vconstraint

k_colors_plot,'main.'+version+'.sav',psfile='main_colors_plot.ps', $
  zrange=[0.,0.23], to_z=0.1,lumlim=[-21.5,-21.2], $
  colorlimits=[[0.5,2.7],[0.23,1.13],[0.25,0.5],[0.01,0.52]], $
  subsample=2, zsplit=[0.05,0.1,0.17], addgrgap=addgrgap, $
  sdssfix=sdssfix, vconstraint=vconstraint, version=version, vpath=vpath
k_colors_plot,'main.'+version+'.sav',psfile='main_colors_plot.z.ps', $
  zrange=[0.,0.23], lumlim=[-21.5,-21.2], $
  colorlimits=[[0.5,2.9],[0.23,1.39],[0.25,0.55],[0.01,0.52]], $
  subsample=2, zsplit=[0.05,0.1,0.17], addgrgap=addgrgap, $
  sdssfix=sdssfix, vconstraint=vconstraint, version=version, vpath=vpath

k_colors_plot,'lrg.'+version+'.sav',psfile='lrg_colors_plot.ps', $
  zrange=[0.1,0.49], to_z=0.3,lumlim=[-22.8,-22.5], $
  nprimtargetmask=67108864, $
  colorlimits=[[0.7,3.3],[0.84,1.74],[0.43,0.68],[0.09,0.59]], $
  zsplit=[0.2,0.3,0.4], addgrgap=addgrgap, $
  sdssfix=sdssfix, vconstraint=vconstraint, version=version, vpath=vpath
k_colors_plot,'lrg.'+version+'.sav',psfile='lrg_colors_plot.z.ps', $
  zrange=[0.1,0.49], lumlim=[-22.8,-22.5], $
  nprimtargetmask=67108864, $
  colorlimits=[[0.5,3.6],[0.84,2.04],[0.39,0.88],[0.09,0.59]], $
  zsplit=[0.2,0.3,0.4], addgrgap=addgrgap, $
  sdssfix=sdssfix, vconstraint=vconstraint, version=version, vpath=vpath

if(not keyword_set(dots)) then dots=[0.08,0.0,-0.15]
k_coeffdist_plot,'all.'+version+'.sav',psfile='k_coeffdist_plot.ps', $
  subsample=3,basecoeff=3,version=version,vpath=vpath, $
  dots=dots, scalespec=[1.0,1.4,1.5],addgrgap=addgrgap, $
  sdssfix=sdssfix, vconstraint=vconstraint
k_kcorrect_plot,'all.'+version+'.sav',psfile='k_kcorrect_plot.ps', $
  subsample=4, to_z=0.3, addgrgap=addgrgap, $
  sdssfix=sdssfix, vconstraint=vconstraint, version=version, vpath=vpath
k_model_plot,'all.'+version+'.sav',version=version,vpath=vpath, $
  psfile='k_model_plot.ps',subsample=3,sdssfix=sdssfix, $
  vconstraint=vconstraint, addgrgap=addgrgap, nsig=15

;setenv,'SPECTRO_DATA=/global/data/sdss/spectro'
;k_speck,'main.'+version+'.sav','all.'+version+'.specK.0.1.sav',subsample=5, to_z=0.1
;$touch gotpast
;k_speck,'lrg.'+version+'.sav','all.'+version+'.specK.0.3.sav',subsample=5, to_z=0.3
k_speck_plot,'all.'+version+'.specK.0.1.sav',to_z=0.1,ylimits=[-0.16,0.16], $
  psfile='k_speck_plot.0.1.ps',zrange=[0.02,0.21],primtargetmask=64, $
  lumlim=[-23.,-19.], version=version, vpath=vpath, $
  addgrgap=addgrgap,vconstraint=vconstraint,sdssfix=sdssfix
k_speck_plot,'all.'+version+'.specK.0.3.sav',to_z=0.3,ylimits=[-0.29,0.29], $
  psfile='k_speck_plot.0.3.ps',zrange=[0.1,0.5], version=version, vpath=vpath
k_speck_plot,'all.'+version+'.specK.0.1.sav',to_z=0.1,ylimits=[-0.16,0.16], $
  psfile='k_speck_plot.fitfib.0.1.ps',zrange=[0.02,0.21],primtargetmask=64, $
  lumlim=[-23.,-19.],/fitfib, version=version, vpath=vpath, $
  addgrgap=addgrgap,vconstraint=vconstraint,sdssfix=sdssfix
k_speck_plot,'all.'+version+'.specK.0.3.sav',to_z=0.3,ylimits=[-0.29,0.29], $
  psfile='k_speck_plot.fitfib.0.3.ps',zrange=[0.1,0.5],/fitfib, $
  version=version, vpath=vpath, $
  addgrgap=addgrgap,vconstraint=vconstraint,sdssfix=sdssfix
k_speck_plot,'all.'+version+'.specK.0.1.sav',to_z=0.1,nsig=3.4, $
  psfile='k_speck_plot.fiber.0.1.ps',usefiber='main.'+version+'.sav', $
  zrange=[0.02,0.21], version=version, vpath=vpath, $
  addgrgap=addgrgap,vconstraint=vconstraint,sdssfix=sdssfix
k_speck_plot,'all.'+version+'.specK.0.3.sav',to_z=0.3,nsig=1.5, $
  psfile='k_speck_plot.fiber.0.3.ps',usefiber='lrg.'+version+'.sav', $
  zrange=[0.1,0.5],lumlim=[-24.0,-21.5], version=version, vpath=vpath, $
  addgrgap=addgrgap,vconstraint=vconstraint,sdssfix=sdssfix

;$cp -f k_kcorrect_plot.ps main_colors_plot.ps lrg_colors_plot.ps k_coeffdist_plot.ps k_speck_plot.0.1.ps k_speck_plot.0.3.ps k_speck_plot.fitfib.0.1.ps k_speck_plot.fitfib.0.3.ps $KCORRECT_DIR/docs/paper

end
