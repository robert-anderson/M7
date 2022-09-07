from system_test import *

skip_unless('many-body basis function', 'fermion (determinant)')
run(assets=['HF_RDMs/FCIDUMP'])
check_nw()
check_rdms()
