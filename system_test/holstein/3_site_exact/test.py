from system_test import *

skip_unless('many-body basis function', 'fermion-boson (determinant-permanent product)')
run()
check_nw()
