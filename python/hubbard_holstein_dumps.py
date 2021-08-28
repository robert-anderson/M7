import numpy

header = ''' &FCI NORB= {}
 &END
'''

def write_ebdump(v, v_unc=None, fname='EBDUMP'):
    '''
    write the coeffients of boson "excitations" and "de-excitations" which correspond to
    single fermion number-conserving excitations (ranksigs 1110, and 1101 in M7 nomenclature)
    and those which do not couple to fermion degrees of freedom at all (0010, and 0001),
    
    these coefficients are given here by the v, and v_unc (uncoupled) args respectively
    '''
    assert len(v.shape)==3
    nsite = v.shape[0]
    if v_unc is None:
        v_unc = numpy.zeros(nsite)
    elif len(numpy.shape(v_unc))==0:
        v_unc = numpy.ones(nsite)*v_unc
    else: assert len(numpy.shape(v_unc))==1
        
    with open(fname, 'w') as f:
        f.write(header.format(nsite))
        for n in range(nsite):
            for p in range(nsite):
                for q in range(nsite):
                    if (v[n, p, q]==0.0): continue
                    f.write('{}    {}    {}    {}\n'.format(v[n, p, q], n+1, p+1, q+1))
        for n in range(nsite):
            if (v_unc[n]==0.0): continue
            f.write('{}    {}    {}    {}\n'.format(v_unc[n], n+1, 0, 0))

def write_bosdump(w, nsite, fname='BOSDUMP'):
    '''
    write the coeffients of boson number-conserving operators (ranksig 0011 in M7 nomenclature)
    '''
    nsite = w.shape[0]
    ndim = len(w.shape)
    if ndim==1: w = numpy.diag(w)
    ndim = len(w.shape)
    assert ndim==2
    with open(fname, 'w') as f:
        f.write(header.format(nsite))
        for n in range(nsite):
            for m in range(nsite):
                if (w[n, m]==0.0): continue
                f.write('{}    {}    {}\n'.format(w[n, m], n+1, m+1))


def setup_hh(nsite, t, U, omega0, g, apbc):
    '''
    set up standard Hubbard-Holstein Hamiltonian
    '''
    idx = numpy.arange(nsite-1)
    tmat_hh = numpy.zeros((nsite,)*2)
    tmat_hh[idx+1,idx] = tmat_hh[idx,idx+1] = -t
    tmat_hh[0,nsite-1] = tmat_hh[nsite-1,0] = t if apbc else -t

    Umat_hh = numpy.zeros((nsite,)*4)
    for i in range(nsite):
        Umat_hh[i,i,i,i] = U

    Vmat_hh = numpy.zeros((nsite,)*3)
    for i in range(nsite):
        Vmat_hh[i,i,i] = g

    Omat_hh = numpy.zeros((nsite,)*2)
    for i in range(nsite):
        Omat_hh[i,i] = omega0

    return tmat_hh, Umat_hh, Vmat_hh, Omat_hh


def delta_zpm(g, nelec, omega0, nsite):
    '''
    get the constant shift in energy owing to the ZPM transformation
    '''
    return (g**2 * nelec**2) / (omega0*nsite)

def make_zpm_coeffs(tmat, nsite, nelec, g, omega0):
    '''
    effect the ZPM transformation on the one-electron integrals, and also return the "uncoupled" boson
    excitation/de-excitation coefficents, and the ZPM energy shift
    '''
    nav = nelec/float(nsite)
    tmat_trans = tmat.copy()
    tmat_trans -= numpy.eye(nsite)*(2.0*nav*g**2)/omega0
    return tmat_trans, -numpy.ones(nsite)*nav*g, delta_zpm(g, nelec, omega0, nsite)


def write_all(tmat, Umat, Vmat, Vmat_unc, Omat, ecore, nelec, fname_fcidump='FCIDUMP', fname_ebdump='EBDUMP', fname_bosdump='BOSDUMP'):
    '''
    include the delta_zpm value in the core energy of the FCIDUMP
    '''
    from pyscf import tools, ao2mo
    nsite = tmat.shape[0]
    tools.fcidump.from_integrals(fname_fcidump, tmat, 
        ao2mo.restore(8,Umat,nsite), nsite, nelec, ecore, 0, [1,]*nsite)
    
    write_ebdump(Vmat, Vmat_unc, fname_ebdump)
    write_bosdump(Omat, fname_bosdump)



if __name__ == '__main__':
    from eb_solve import kernel
    '''
    integral file writing example and a rough check that the difference between the standard HH and the 
    transformed hamiltonian vanish in the limit of a large boson cutoff
    '''
    nsite, nelec = 4, 4
    t = 1
    U = 4
    omega0 = 0.5
    g = 0.1
    apbc = True

    tmat, Umat, Vmat, Omat = setup_hh(nsite, t, U, omega0, g, apbc)
    write_all(tmat, Umat, Vmat, None, Omat, 0.0, nelec, 'FCIDUMP', 'EBDUMP', 'BOSDUMP')

    tmat_trans, Vmat_unc, delta_zpm = make_zpm_coeffs(tmat, nsite, nelec, g, omega0)
    write_all(tmat_trans, Umat, Vmat, Vmat_unc, Omat, delta_zpm, nelec, 'FCIDUMP_ZPM', 'EBDUMP_ZPM', 'BOSDUMP_ZPM')


    errs = []
    for nboson_max in range(7):
        res0 = kernel(tmat, Umat, Vmat, Omat, nsite, nelec, nsite, nboson_max, adj_zero_pho = False)
        e1 = res0[0]
        print(f'*\n* nboson_max = {nboson_max}, standard energy = {e1}\n*')
        res0 = kernel(tmat_trans, Umat, Vmat, Omat, nsite, nelec, nsite, nboson_max, adj_zero_pho = True)
        e2 = res0[0]+delta_zpm
        print(f'*\n* nboson_max = {nboson_max}, ZPM-removed energy = {e2}\n*')
        errs.append(abs(e1-e2))

    print(errs)
    assert numpy.all(numpy.diff(errs) < 0)
    assert errs[-1]<1e-8

