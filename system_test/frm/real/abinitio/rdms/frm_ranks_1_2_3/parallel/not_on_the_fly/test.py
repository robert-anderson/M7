from system_test import *

require_mbf_type('fermion')
require_ham_arith('real')

run(link_deps=['HF_RDMs/FCIDUMP'], nrank=2)

compare_nw()
compare_rdm_archives()

exact_e0 = -99.94213890
check_shift(exact_e0)
check_proje(exact_e0)

check_spinfree_rdms('M7.rdm.h5', 'HF_RDMs/exact_rdms.pkl', ('1100', '2200', '3300'))
