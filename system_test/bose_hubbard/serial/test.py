from system_test import *

skip_unless('many-body basis function', 'boson (permanent)')
run()
check_nw()
exact_e0 = -4.40663255
check_shift_correct(exact_e0, 1000)
check_proje_correct(exact_e0, 1000)
