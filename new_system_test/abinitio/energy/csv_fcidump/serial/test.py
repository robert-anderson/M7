from system_test import *

#skip_if('many-body basis function', 'boson (permanent)')
run(link_deps=['HF_RDMs/FCIDUMP'])
#exact_e0 =-99.9421389039332
comp_nw()
#check_ninit()
#check_shift_correct(exact_e0, 10000)
#check_proje_correct(exact_e0, 10000)
