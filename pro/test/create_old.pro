pro create_old, maggies, maggies_ivar, redshift, kcorrect_old, $
                band_shift=band_shift

if(NOT keyword_set(band_shift)) then band_shift=0.

spawn,'pwd',cwd
spawn,'mkdir -p tmp'
cd,'tmp'
;spawn,'cvs -d blanton@spectro.princeton.edu:/usr/local/cvsroot '+ $
;  'co -r v1_17 kcorrect'
cd,'kcorrect'
spawn,'setup -r '+cwd+'/tmp/kcorrect kcorrect ; kevilmake -k all'
print,'setup -r '+cwd+'/tmp/kcorrect kcorrect'
cd,'../../'

save,maggies,maggies_ivar,redshift,filename='tmp.sav'
spawn,'setup  -r '+cwd+'/tmp/kcorrect kcorrect ; ' + $
  ' echo restore,\"tmp.sav\" \& '+ $
  ' kcorrect, maggies, maggies_ivar, redshift, kcorrect_old, /invvar, /maggies, kcorrectz='+strtrim(string(band_shift),2)+' \& '+ $
  ' save,filename=\"tmp_kcorrect_old.sav\" | idl'

restore,'tmp_kcorrect_old.sav'

end

restore,'testdata*.sav'
test_vs_old,maggies[0:4,0:n_elements(gals)-1L], $
  maggies_ivar[0:4,0:n_elements(gals)-1L], $
  gals.redshift
  
