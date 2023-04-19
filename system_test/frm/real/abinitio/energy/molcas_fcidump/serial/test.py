from system_test import *

require_mbf_type('fermion')
require_ham_arith('real')

run(link_deps=[('N2_Molcas/molcas.FciDmp.h5', 'FCIDUMP')])

exact_e0 = -109.02180323
compare_nw()
check_shift(exact_e0)
check_proje(exact_e0)
