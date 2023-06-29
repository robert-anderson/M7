from system_test import *

require_mbf_type('fermion')
require_ham_arith('real')

run('config.large_ci.yaml', link_deps=['HF_RDMs/FCIDUMP'], nrank=2)
run('config.rdms.yaml', link_deps=['HF_RDMs/FCIDUMP'], nrank=2)

exact_e0 = -99.94213890

check_rdm_energy(exact_e0)
