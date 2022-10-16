from system_test import *

skip_unless('many-body basis function', 'fermion (determinant)')
run(assets=[('N2_Molcas/molcas.FciDmp.h5', 'FCIDUMP'), ('N2_Molcas/f_act.h5', 'fock.h5')])
check_nw()
check_rdm_archives()
