from system_test import *

skip_if('many-body basis function', 'boson (permanent)')
run(nrank=1, assets=['N_UHF/FCIDUMP.spin_blocks'])
exact_e0 =-54.419939662413
check_nw()
check_ninit()
check_shift_correct(exact_e0, 30000)
check_proje_correct(exact_e0, 30000)
