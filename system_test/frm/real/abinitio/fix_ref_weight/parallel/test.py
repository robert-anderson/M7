from system_test import *

require_mbf_type('fermion')
require_ham_arith('real')

run(nrank=2, link_deps=['RHF_N2_6o6e/FCIDUMP'])
compare_nw()
compare_shift()
compare_ninit()
compare_nocc_mbf()
