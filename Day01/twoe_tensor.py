
# here define the initial tensor for the two-electron repulsion integrals (ERI)
Vee = zeros((nel,nel,nel,nel))
# this is a very specific way to generate the ERI and is not generalized
for i in range(0, N):
    for j in range(0, N):
        for k in range(0, N):
            for l in range(0, N):
                # compute R distances for each center
                RAP = A2[i]*R/(A2[i] + A1[j])
                RBP = R - RAP
                RAQ = A2[k]*R/(A2[k] + A1[l])
                RBQ = R - RAQ
                RPQ = RAP - RAQ
                Vee[0,0,0,0] += twoe(A1[i], A1[j], A1[k], A1[l], 0.0, 0.0,    0.0)*D1[i]*D1[j]*D1[k]*D1[l]
                Vee[1,0,0,0] += twoe(A2[i], A1[j], A1[k], A1[l],  R2, 0.0, RAP**2)*D2[i]*D1[j]*D1[k]*D1[l]
                Vee[1,0,1,0] += twoe(A2[i], A1[j], A2[k], A1[l],  R2,  R2, RPQ**2)*D2[i]*D1[j]*D2[k]*D1[l]
                Vee[1,1,0,0] += twoe(A2[i], A2[j], A1[k], A1[l], 0.0, 0.0,     R2)*D2[i]*D2[j]*D1[k]*D1[l]
                Vee[1,1,1,0] += twoe(A2[i], A2[j], A2[k], A1[l], 0.0,  R2, RBQ**2)*D2[i]*D2[j]*D2[k]*D1[l]
                Vee[1,1,1,1] += twoe(A2[i], A2[j], A2[k], A2[l], 0.0, 0.0,    0.0)*D2[i]*D2[j]*D2[k]*D2[l]

# Fill the rest of the ERI using symmetry relations 
Vee[0,1,0,0] = Vee[1,0,0,0] 
Vee[0,0,1,0] = Vee[1,0,0,0] 
Vee[0,0,0,1] = Vee[1,0,0,0]
Vee[0,1,1,0] = Vee[1,0,1,0]
Vee[1,0,0,1] = Vee[1,0,1,0]
Vee[0,1,0,1] = Vee[1,0,1,0]
Vee[0,0,1,1] = Vee[1,1,0,0]
Vee[1,1,0,1] = Vee[1,1,1,0]
Vee[1,0,1,1] = Vee[1,1,1,0]
Vee[0,1,1,1] = Vee[1,1,1,0]
