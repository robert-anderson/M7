from hubbard_holstein_dumps import *
from eb_solve import kernel

#w=0.5, U=4, g=sqrt(0.15)
nsite, nelec = 8, 8
t = 1
U = 4
omega0 = 0.5
g = 0.15**0.5
apbc = True
nboson_max = 1

tmat, Umat, Vmat, Omat = setup_hh(nsite, t, U, omega0, g, apbc)
write_all(tmat, Umat, Vmat, None, Omat, 0.0, nelec, 'FCIDUMP', 'EBDUMP', 'BOSDUMP')

tmat_trans, Vmat_unc, delta_zpm = make_zpm_coeffs(tmat, nsite, nelec, g, omega0)
write_all(tmat_trans, Umat, Vmat, Vmat_unc, Omat, delta_zpm, nelec, 'FCIDUMP_ZPM', 'EBDUMP_ZPM', 'BOSDUMP_ZPM')

res0 = kernel(tmat_trans, Umat, Vmat, Omat, nsite, nelec, nsite, nboson_max, adj_zero_pho = True)
e2 = res0[0]+delta_zpm
print(e2)
#print(f'*\n* nboson_max = {nboson_max}, ZPM-removed energy = {e2}\n*')
