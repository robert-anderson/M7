from system_test import *

skip_if('many-body basis function', 'boson (permanent)')
run(nrank=2, assets=['RHF_N2_6o6e/FCIDUMP'])
check_nw()
check_ninit()
check_nocc_mbf()
