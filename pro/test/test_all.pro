; test all features
pro test_all

failed=0
test_kcorrect,failed
if(failed) then klog,'really bad failure of test_kcorrect'

failed=0
test_kcorrect_standalone,failed
if(failed) then klog,'really bad failure of test_kcorrect_standalone'

failed=0
test_kphotoz,failed
if(failed) then klog,'really bad failure of test_kphotoz'

failed=0
test_kphotoz_standalone,failed
if(failed) then klog,'really bad failure of test_kphotoz_standalone'

end
