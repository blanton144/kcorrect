function lf_test_root_func, px

return,(px-3.)^3

end

pro lf_test_root

print,zbrent(1.,2.,func_name='lf_test_root_func')

end
