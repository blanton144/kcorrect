;; formulae from Bullock et al. (1999) astro-ph/9908159
;; only implemented for z=0; input mass is mass/h
function m2vmax, in_mass

littleh=0.7
omega0=0.3
mass=in_mass/littleh
mstar=1.5e+13/littleh

xx=omega0-1.
delvir=(18.*!DPI^2+82.*xx-39.*xx^2)/omega0
rhou=omega0*(2.776e+11)*littleh^2 ;; in Msolar Mpc^{-3}
rvir=(mass*3./(4.*!DPI)/delvir/rhou)^(1./3.)  ;; in Mpc

gnewton_si=6.67e-11 ;; in m^3 kg^{-1} s^{-2}
kmperm=0.001
mpcperm=1./(3.086e+22)
kgpermsolar=1.989e+30
gnewton=gnewton_si*(kmperm)^2*mpcperm*kgpermsolar 
;; in Mpc (m s^{-1})^2 Msolar^{-1} 

vvir=sqrt(gnewton*mass/rvir)

mu=mass/mstar
cvir=9.*mu^(-0.13)

acvir=alog(1.+cvir)-cvir/(1.+cvir)

vmax=vvir*sqrt(0.216*cvir/acvir)

return,vmax

end
