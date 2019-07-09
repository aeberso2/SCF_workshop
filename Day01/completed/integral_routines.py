
# the overlap  integral for unnormalized primitives
def S(A, B, RAB2):
    
    S_int = (pi/(A + B))**1.5*exp(-A*B*RAB2/(A + B))
    
    return S_int

# the kinetic energy integrals for unnormalized primitives
def T(A, B, RAB2):
    
    T_int = A*B/(A + B)*(3 - 2*A*B*RAB2/(A + B))*(pi/(A + B))**1.5*exp(-A*B*RAB2/(A + B))
    
    return T_int

# the potential energy nuclear attraction integrals, 
# we use F0 here and this will only hold for s-type functions
def V(A, B, RAB2, RCP2, ZC):
    
    V_int = 2*pi/(A + B)*F0((A + B)*RCP2)*exp(-A*B*RAB2/(A + B))
    V_int *= -ZC
    
    return V_int

# the two electron integral for unnormalized primitives 
# A, B, C, D are the exponents alpha, beta, etc. 
# RAB2 = squared distance bewteen center A and center B, etc. 
# Again, F0 is used so DON'T USE IT FOR P-orbitals!!
def twoe(A, B, C, D, RAB2, RCD2, RPQ2):
    
    twoe_int = (2*pi**2.5)/((A + B)*(C + D)*sqrt(A + B + C + D))
    twoe_int *= F0((A + B)*(C + D)*RPQ2/(A + B + C + D))
    twoe_int *= exp(-A*B*RAB2/(A + B) - C*D*RCD2/(C + D))
    return twoe_int

# form the G-matrix from the electron-electron repulsion and density matrices
def formG(P):
    g = zeros((nel,nel))
    for i in range(nel):
        for j in range(nel):
            for k in range(nel):
                for l in range(nel):
                    # Coulomb integral -> J
                    g[i,j] += 2*P[k,l]*Vee[i,j,k,l]
                    # exchange integral -> K
                    g[i,j] -= P[k,l]*Vee[i,k,j,l]
    return g
