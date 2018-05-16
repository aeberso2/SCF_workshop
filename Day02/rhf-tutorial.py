import psi4
from numpy import *
import pandas as pd
import os, sys
from scipy.constants import physical_constants
from scipy.linalg import eigh, inv

bohr2ang = physical_constants['Bohr radius'][0]*1e10
os.system('rm -f output.dat timer.dat')

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
    B[-1, :] = -1
    B[ :,-1] = -1
    B[-1,-1] = 0
    for i in range(len(F_list)):
        for j in range(len(F_list)):
            B[i,j] = sum(diis_resid[i]*diis_resid[j])

    rhs = zeros((B_dim))
    rhs[-1] = -1

    coeff = linalg.solve(B, rhs)

    F_diis = zeros_like(F_list[0])
    for ix in range(coeff.shape[0] - 1):
        F_diis += coeff[ix]*F_list[ix]

    return F_diis

def formG(P):
    J = einsum('pqrs,rs->pq', Vee, P, optimize=True)
    K = einsum('prqs,rs->pq', Vee, P, optimize=True)
    G = 2*J - K
    return G

# psi4 introduction!

# first we can set an outfile if we want to dump all the calc info
psi4.core.set_output_file('output.dat')

bond_dist = 1.4632*bohr2ang

cmpd = 'HeH'

# here is how we define the geometry
mol = psi4.geometry("""
        He
        H 1 {: .5f}
        symmetry c1
        """.format(bond_dist))

# set charge to positive 1
mol.set_molecular_charge(1)

the_basis = 'sto-3g'
# set our calculation options
psi4.set_options({'guess':'core',
                  'basis':'{}'.format(the_basis),
                  'scf_type':'pk',
                  'e_convergence':1e-8,
                  'reference': 'rhf'})

# our guess is the core guess, -> can also be sad (superposition of atomic densities)
# scf_type -> ERI algorithm, pk is default
# reference -> rhf, uhf, and maybe rohf

# compute static 1e- and 2e- quantities in Psi4
# Class initialization
wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option('basis'))
# mints is the integral helper
mints = psi4.core.MintsHelper(wfn.basisset())

# the Smat is the atomic orbital overlap
Smat = asarray(mints.ao_overlap())
# number of basis functions, alpha orbitals -> rhf so just call alpha
nbf = wfn.nso()
ndocc = wfn.nalpha()

# Build core Hamiltonian
Tmat = asarray(mints.ao_kinetic())
Vmat = asarray(mints.ao_potential())
Hmat = Tmat + Vmat

# build the nasty two-electron repulsion integral
Vee = asarray(mints.ao_eri())

# Construct AO orthogonalization matrix A
# this is the Psi4 way, which is for symmetric orthog
# A = mints.ao_overlap()
# A.power(-0.5, 1.0e-16)
# A = asarray(A)

# get nuclear repulsion energy from Psi4
E_nuc = mol.nuclear_repulsion_energy()


# we'll keep our way in here, it works the same
u, V = eigh(Smat)
U = sqrt(inv(u*eye(len(u))))
A = dot(V.T, dot(U, V))


# maximum scf iterations
maxiter = 40

# energy convergence criterion
E_conv = 1.0e-6
D_conv = 1.0e-4

# pre-iteration step
# scf & previous energy
SCF_E = 0.0
E_old = 0.0
# form core guess
C, P, epsilon = diag_F(Hmat, ndocc)

# trial and resiual vector lists
F_list = []
R_list = []
diis_resid = []

print('Number of occupied orbitals {}'.format(ndocc))
print('Number of basis functions {}'.format(nbf))

print('==> Starting SCF Iterations <==\n')

# comment in to write initial wavefunction
# f = open('init_wfn.dat', 'w')
# for i in range(C.shape[0]):
#     for j in range(C.shape[1]):
#         print('{: 23.15f}'.format(C[i,j]), file=f)
# f.close()

for scf_iter in range(maxiter):
    # Build the Fock matrix
    # We will build the G matrix in a slightly different way, using
    # the einsum function
    F = Hmat + formG(P)

    # for the diis
    # A * (F*P*S - S*P*F) * A
    M = dot(F, dot(P, Smat)) - dot(Smat, dot(P, F))
    diis_r = dot(A, dot(M, A))

    F_list.append(F)
    R_list.append(diis_r)

    SCF_E = sum((Hmat + F)*P) + E_nuc

    dE = SCF_E - E_old

    dRMS = mean(diis_r**2)**0.5
    print('SCF Iteration {:3d}: Energy = {: 4.16f} dE = {: 1.5e} dRMS = {:1.5e}'.format(scf_iter+1, SCF_E, dE, dRMS))

    if (abs(dE) < E_conv) and (dRMS < D_conv):
        break
    E_old = SCF_E

    if scf_iter >= 2:
        F = diis_xtrap(F_list, R_list)

    C, P, epsilon = diag_F(F, ndocc)

    if scf_iter == maxiter:
        psi4.core.clean()
        raise Exception("Maximum number of SCF iterations exceeded.")

print('\nSCF Converged.')
print('Final RHF Energy: {: .8f} [Eh]'.format(SCF_E))

# print the final wavefunction
f = open('final_wfn.dat', 'w')
for i in range(C.shape[0]):
    for j in range(C.shape[1]):
        print('{: 23.15f}'.format(C[i,j]), file=f)
f.close()

SCF_E_psi, wfn = psi4.energy('SCF', return_wfn=True)
psi4.compare_values(SCF_E_psi, SCF_E, 6, 'SCF Energy')
# remove any molden file with our compoudns name if it exists
os.system('rm -f {}.rhf.molden'.format(cmpd))
# create new molden file
psi4.molden(wfn, '{}.rhf.molden'.format(cmpd))

# uncomment for FCI energy
# FCI_E_psi = psi4.energy('FCI')
# print(FCI_E_psi)
