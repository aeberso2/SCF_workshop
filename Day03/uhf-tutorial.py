import psi4
from numpy import *
import pandas as pd
import os, sys
from numpy.linalg import eigh, inv
from scipy.constants import physical_constants
from scipy.linalg import eigh, inv

bohr2ang = physical_constants['Bohr radius'][0]*1e10

os.system('rm -f output.dat timer.dat')

# these can stay the same
def diag_F(F, norb):
    Fp = dot(A, dot(F, A))
    e, Cp = linalg.eigh(Fp)
    C = dot(A, Cp)
    C_occ = C[:, :norb]
    P = einsum('pi,qi->pq', C_occ, C_occ)
    return (C, P, e)


def diis_xtrap(F_list, diis_resid):
    B_dim = len(F_list) + 1
    B = empty((B_dim, B_dim))
    B[-1,  :] = -1
    B[ :, -1] = -1
    B[-1, -1] = 0
    for i in range(len(F_list)):
        for j in range(len(F_list)):
            B[i,j] = einsum('ij,ij->', diis_resid[i], diis_resid[j])

    rhs = zeros((B_dim))
    rhs[-1] = -1

    coeff = linalg.solve(B, rhs)
    F_diis = zeros_like(F_list[0])
    for ix in range(coeff.shape[0] - 1):
        F_diis += coeff[ix]*F_list[ix]

    return F_diis


psi4.core.set_output_file('output.dat', False)

bond_dist = 1.4632*bohr2ang

cmpd = 'HeH'

# here is how we define the geometry
mol = psi4.geometry("""
        He
        H 1 {: .5f}
        symmetry c1
        """.format(bond_dist))

the_basis = 'sto-3g'
psi4.set_options({'guess':'core',
                  'basis':'{}'.format(the_basis),
                  'scf_type':'pk',
                  'e_convergence':1e-8,
                  'reference': 'uhf'})

mol.set_molecular_charge(1)

maxiter = 40
# energy convergence criterion
E_conv = 1.0e-6
D_conv = 1.0e-4

# compute static 1e- and 2e- quantities in Psi4
wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option('basis'))
mints = psi4.core.MintsHelper(wfn.basisset())

Smat = asarray(mints.ao_overlap())
# number of basis functions, alpha & beta orbitals and doubly occupied orbitals
nbf = wfn.nso()
nalpha = wfn.nalpha()
nbeta = wfn.nbeta()
ndocc = min(nalpha, nbeta)

# build the nasty two-eri tensor again...
Vee = asarray(mints.ao_eri())

# Build core hamiltonian
Tmat = asarray(mints.ao_kinetic())
Vmat = asarray(mints.ao_potential())
Hmat = Tmat + Vmat

# Construct AO orthogonalization matrix A
# A = mints.ao_overlap()
# A.power(-0.5, 1.0e-16)
# A = asarray(A)

# orthogonalize mit symmetric orthog
u, V = eigh(Smat)
U = sqrt(inv(u*eye(len(u))))
A = dot(V, dot(U, V.T))

# build alpha & beta core guess
# must do it for two elements
Ca, Pa, epa =
Cb, Pb, epb =

# get nuclear repulsion energy
E_nuc = mol.nuclear_repulsion_energy()

# pre-iteration step
# scf & previous energy
SCF_E = 0.0
E_old = 0.0

# trial and resiual vector lists -- one for each alpha + beta
F_list_a = []
F_list_b = []
R_list_a = []
R_list_b = []

print('Number of basis functions {}'.format(nbf))
print('Number of singly occupied orbitals {}'.format(abs(nalpha - nbeta)))
print('Number of doubly occupied orbitals {}'.format(ndocc))

print('==> Starting SCF Iterations <==\n')

# everything will be in twos
for scf_iter in range(maxiter):

    # Build Fa and Fb matrices
    # we have make the Ja and Jb forms
    Ja =
    Jb =

    Ka =
    Kb =

    Fa =
    Fb =

    # build two M's this time
    Ma =
    Mb =
    diis_r_a =
    diis_r_b =

    F_list_a.append(Fa)
    F_list_b.append(Fb)

    R_list_a.append(diis_r_a)
    R_list_b.append(diis_r_b)

    SCF_E =

    dE = SCF_E - E_old
    dRMS = 0.5*(mean(diis_r_a**2)**0.5 + mean(diis_r_b**2)*0.5)

    print('SCF Iteration {:3d}: Energy = {: 4.16f} dE = {: 1.5e} dRMS = {:1.5e}'.format(scf_iter + 1, SCF_E, dE, dRMS))


    if (abs(dE) < E_conv) and (dRMS < D_conv):
        break

    E_old = SCF_E

    if scf_iter >= 2:
        Fa = diis_xtrap(F_list_a, R_list_a)
        Fb = diis_xtrap(F_list_b, R_list_b)

    Ca, Pa, epa =
    Cb, Pb, epb =

    if scf_iter == maxiter:
        psi4.core.clean()
        raise Exception("Maximum number of SCF iterations exceeded.")

print('\nSCF Converged.')
print('Final UHF Energy: {: .8f} [Eh]'.format(SCF_E))

SCF_E_psi, wfn = psi4.energy('SCF', return_wfn=True)
psi4.compare_values(SCF_E_psi, SCF_E, 6, 'SCF Energy')

# create molden file
os.system('rm -f {}.uhf.molden'.format(cmpd))
psi4.molden(wfn, '{}.uhf.molden'.format(cmpd))
