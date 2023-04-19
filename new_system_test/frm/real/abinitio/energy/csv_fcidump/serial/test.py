from system_test import *

require_mbf_type('fermion')
require_ham_arith('real')

run(link_deps=['HF_RDMs/FCIDUMP'])
exact_e0 =-99.9421389039332

compare_nw()
compare_ninit()

check_shift(exact_e0, 10000)
check_proje(exact_e0, 10000)
