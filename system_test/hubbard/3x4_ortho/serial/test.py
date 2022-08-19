from system_test import *

skip_unless('many-body basis function', 'fermion (determinant)')
run()
check_nw()
