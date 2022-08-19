from system_test import *

skip_unless('many-body basis function', 'fermion (determinant)')
run(nrank=2)
check_nw()
