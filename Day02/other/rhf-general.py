import psi4
from numpy import *
import pandas as pd
import os, sys
from scipy.constants import physical_constants
import matplotlib.pyplot as plt
from scipy.linalg import eigh, inv

bohr2ang = physical_constants['Bohr radius'][0]*1e10
os.system('rm -f output.dat timer.dat')


cmpd = 'HeH'
# cmpd = 'H2'
# psi4 introduction!
def run_scf(bond_dist, out_file):

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

    psi4.core.set_output_file('output.dat', False)

    mol = psi4.geometry("""
            He
            H 1 {: .5f}
            symmetry c1
            """.format(bond_dist))
    # set charge to positive 1
    mol.set_molecular_charge(1)
    mol.set_multiplicity(1)
    the_basis = 'sto-3g'
#     the_basis = '3-21G'
#     the_basis = '6-31G*'
#     the_basis = 'cc-pvdz'
    # set our calculation options
    psi4.set_options({'guess':'core',
                      'basis':'{}'.format(the_basis),
                      'scf_type':'pk',
                      'e_convergence':1e-8,
                      'reference': 'rhf'})

    wfn = psi4.core.Wavefunction.build(mol, psi4.core.get_global_option('basis'))
    mints = psi4.core.MintsHelper(wfn.basisset())
    Smat = asarray(mints.ao_overlap())
    nbf = wfn.nso()
    ndocc = wfn.nalpha()

    Tmat = asarray(mints.ao_kinetic())
    Vmat = asarray(mints.ao_potential())
    Hmat = Tmat + Vmat

    Vee = asarray(mints.ao_eri())
    E_nuc = mol.nuclear_repulsion_energy()
#     u, V = eigh(Smat)
#     U = sqrt(inv(u*eye(len(u))))
#     A = dot(V.T, dot(U, V))

#     u, V = eigh(Smat)
#     u = 1/sqrt(u)
#     A = V*u

    u, V = eigh(Smat)
    U = sqrt(inv(u*eye(len(u))))
    A = dot(V, dot(U, V.T))
    AT = A.T

    maxiter = 40
    E_conv = 1.0e-6
    D_conv = 1.0e-4

    SCF_E = 0.0
    E_old = 0.0
    C, P, epsilon = diag_F(Hmat, ndocc)

    F_list = []
    R_list = []
    diis_resid = []

    for scf_iter in range(maxiter):
        F = Hmat + formG(P)

        M = dot(F, dot(P, Smat)) - dot(Smat, dot(P, F))
        diis_r = dot(AT, dot(M, A))

        F_list.append(F)
        R_list.append(diis_r)

        SCF_E = sum((Hmat + F)*P) + E_nuc

        dE = SCF_E - E_old

        dRMS = mean(diis_r**2)**0.5

        if (abs(dE) < E_conv) and (dRMS < D_conv):
            break
        E_old = SCF_E

        if scf_iter >= 2:
            F = diis_xtrap(F_list, R_list)

        C, P, epsilon = diag_F(F, ndocc)

        if scf_iter == maxiter:
            psi4.core.clean()
            raise Exception("Maximum number of SCF iterations exceeded.")

    print('{: .5f}  {: .8f}'.format(SCF_E, bond_dist), file=out_file)
    return SCF_E





# SCF_E_psi, wfn = psi4.energy('SCF', return_wfn=True)
# psi4.compare_values(SCF_E_psi, SCF_E, 6, 'SCF Energy')
# remove any molden file with our compoudns name if it exists
os.system('rm -f {}.rhf.molden'.format(cmpd))
# create new molden file
# psi4.molden(wfn, '{}.rhf.molden'.format(cmpd))

bond_dist = 1.4632*bohr2ang
outfile = open('pes.HeH+.dat', 'w')
r = linspace(1.05, 7.5, 60)*bohr2ang
pesvals = zeros(60)
for i in range(len(r)):
    pesvals[i] = run_scf(r[i], outfile)

fig, ax = plt.subplots()
ax.plot(r/bohr2ang, pesvals)
plt.show()




# uncomment for FCI energy
# FCI_E_psi = psi4.energy('FCI')
# print(FCI_E_psi)
