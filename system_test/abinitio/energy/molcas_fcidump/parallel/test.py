from system_test import *

skip_if('many-body basis function', 'boson (permanent)')
run(nrank=2, assets=[('N2_Molcas/molcas.FciDmp.h5', 'FCIDUMP')])
exact_e0 = -109.02180323
check_nw()
check_shift_correct(exact_e0, 12500)
check_proje_correct(exact_e0, 12500)
