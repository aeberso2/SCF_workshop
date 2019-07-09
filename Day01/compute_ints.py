
# Calculate the one-electron integrals: -> overlap, kinetic, and potential 
# center A is the first atom, center B the second -> origin is on A
# sum over all N-indices for each atom 
for i in range(0, N):
    for j in range(0, N):
        RAP = A2[j]*R/(A1[i] + A2[j])
        # RAP2 is squared distance between center A and center P, etc. 
        RAP2 = RAP**2
        RBP2 = (R - RAP)**2
        Smat[0,1]  += S(A1[i], A2[j], R2)*D1[i]*D2[j]
        # note the symmetry, the lower off-diagonal elements are equal to the upper off-diagonal elements
        Smat[1,0] = Smat[0,1]
        
        Tmat[0,0] += T(A1[i], A1[j], 0.0)*D1[i]*D1[j]
        Tmat[0,1] += T(A1[i], A2[j],  R2)*D1[i]*D2[j]
        Tmat[1,0] = Tmat[0,1]
        Tmat[1,1] += T(A2[i], A2[j], 0.0)*D2[i]*D2[j]
        
        # Potential from atom A + potential from atom B 
        Vmat[0,0] += V(A1[i], A1[j], 0.0, 0.0, ZA)*D1[i]*D1[j]
        Vmat[0,0] += V(A1[i], A1[j], 0.0,  R2, ZB)*D1[i]*D1[j]
        
        # off-diagonal nuclear attraction to center A, etc.
        Vmat[0,1] += V(A1[i], A2[j], R2, RAP2, ZA)*D1[i]*D2[j] 
        Vmat[0,1] += V(A1[i], A2[j], R2, RBP2, ZB)*D1[i]*D2[j]

        # note the symmetry, the off-diagonal elements all have the same value 
        Vmat[1,0] = Vmat[0,1]
        
        Vmat[1,1] += V(A2[i], A2[j], 0.0,  R2, ZA)*D2[i]*D2[j]
        Vmat[1,1] += V(A2[i], A2[j], 0.0, 0.0, ZB)*D2[i]*D2[j]

# set diagonal of Smat to one
Smat[0,0] = 1.0
Smat[1,1] = 1.0

# create the core hamiltonian, H = T + V
Hmat = Tmat + Vmat 
