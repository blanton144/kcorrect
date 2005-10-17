pro deep_stack

deep=mrdfits(getenv('KCORRECT_DIR')+ $
             '/data/redshifts/deep/zcat.dr1.uniq.fits.gz',1)
igot=where(deep.zquality ge 3 AND $
           deep.zhelio gt 0.01 AND $
           abs(deep.dec) lt 2.)
help,igot
isort=deep[igot].magr
stack_and_measure, deep[igot[isort]].ra, deep[igot[isort]].dec, rerun=137, $
  fl=fl

save, filename='data_deep_stack.sav'

end
