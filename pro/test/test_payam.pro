pro test_payam

readcol, 'mags.txt', comment=';', specz, u, g, r, i, z, $
  uerr, gerr, rerr, ierr, zerr

uerr=sqrt(uerr^2+0.05^2)
gerr=sqrt(gerr^2+0.02^2)
rerr=sqrt(rerr^2+0.02^2)
ierr=sqrt(ierr^2+0.02^2)
zerr=sqrt(zerr^2+0.03^2)

u=u[0:999]
g=g[0:999]
r=r[0:999]
i=i[0:999]
z=z[0:999]
uerr=uerr[0:999]
gerr=gerr[0:999]
rerr=rerr[0:999]
ierr=ierr[0:999]
zerr=zerr[0:999]

mags=transpose(reform([u,g,r,i,z], n_elements(u), 5))
magserr=transpose(reform([uerr,gerr,rerr,ierr,zerr], n_elements(u), 5))
kphotoz, mags, magserr, photoz, /mag, /std

end
