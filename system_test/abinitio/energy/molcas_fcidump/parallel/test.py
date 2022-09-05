from system_test import *

skip_unless('many-body basis function', 'fermion (determinant)')
run(nrank=2, assets=[('N2_Molcas/molcas.FciDmp.h5', 'FCIDUMP')])
exact_e0 = -109.02180323
check_nw()
check_shift_correct(exact_e0, 10000)
check_proje_correct(exact_e0, 10000)
