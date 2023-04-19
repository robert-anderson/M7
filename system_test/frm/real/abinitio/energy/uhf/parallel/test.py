from system_test import *

require_mbf_type('fermion')
require_ham_arith('real')

run(nrank=2, link_deps=['N_UHF/FCIDUMP.spin_blocks'])
exact_e0 =-54.419939662413
compare_nw()
compare_ninit()
check_shift(exact_e0)
check_proje(exact_e0)
