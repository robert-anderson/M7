from system_test import *

skip_if('many-body basis function', 'boson (permanent)')
run('config.accum.yaml', assets=['RHF_N2_6o6e/FCIDUMP'])
run('config.perma.yaml', assets=['RHF_N2_6o6e/FCIDUMP'])
check_nw()
