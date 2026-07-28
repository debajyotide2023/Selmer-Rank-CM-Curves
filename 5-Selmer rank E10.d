#############################################################
# Compute dim Cl(L_d)(chi_d)
# with CM normalization tau(P)=2P
#############################################################

reset()

d = -3


#############################################################
# K = Q(sqrt(-11))
#############################################################

K.<a> = QuadraticField(-19)

R.<X> = PolynomialRing(K)


#############################################################
# F = K(beta)
#############################################################

b1 = -2*a - 38
b2 = (152*a + 1672)/5

F.<beta> = K.extension(X^2+b1*X+b2)


#############################################################
# Elliptic curve
#############################################################

A = -608
B = 5776



#############################################################
# gamma_d^2 = d^3 r(beta)
#############################################################

r = (A+b1^2-b2)*beta + (B+b1*b2)

S.<Y> = PolynomialRing(F)

Ld.<gamma> = F.extension(Y^2-d^3*r)



print("Relative field:")
print(Ld)



#############################################################
# Absolute field + explicit embedding
#############################################################

Labs = Ld.absolute_field('z')

phi = Ld.embeddings(Labs)[0]


a_abs = phi(a)
beta_abs = phi(beta)
gamma_abs = phi(gamma)


print("Check K embedding:")
print(a_abs^2)
print(beta_abs^2 + phi(b1)*beta_abs + phi(b2))


#############################################################
# Check elliptic equation
#############################################################

print("Check gamma relation:")
print(gamma_abs^2 == phi(d^3*r))



#############################################################
# Class group
#############################################################

Cl = Labs.class_group()

print("Class group:")
print(Cl)



#############################################################
# Extract 5-torsion
#############################################################

five_candidates=[]

for g in Cl.gens():

    n=g.order()

    if n % 5 == 0:

        h=g^(n//5)

        if h.order()==5:
            five_candidates.append(h)



def in_span(h,basis):

    r=len(basis)

    for v in GF(5)^r:

        x=Cl(1)

        for i in range(r):
            x *= basis[i]^Integer(v[i])

        if x==h:
            return True

    return False



basis=[]

for g in five_candidates:

    if not in_span(g,basis):
        basis.append(g)



r5=len(basis)


print("5-rank =",r5)



#############################################################
# Define CM automorphism tau(P)=2P
#############################################################

c = (-27*a/76 + 3/4)*beta + 7*a -30


target_beta = phi(-b1-beta)
target_gamma = phi(c*gamma)


Aut = Labs.automorphisms()


tau=None


for s in Aut:

    # tau must lie in Gal(L_d/K)
    # i.e. it fixes sqrt{-11}
    if s(a_abs) == a_abs:

        if s(beta_abs)==target_beta and s(gamma_abs)==target_gamma:
            tau=s
            break



if tau is None:
    raise RuntimeError("CM normalized tau not found")


print("tau found:")
print(tau)


print("\nCheck tau fixes sqrt{-11}:")
print(tau(a_abs) == a_abs)

print("tau(a) =")
print(tau(a_abs))

print("a =")
print(a_abs)



#############################################################
# Coordinate function
#############################################################

def coordinates(h,basis):

    r=len(basis)

    for v in GF(5)^r:

        x=Cl(1)

        for i in range(r):
            x *= basis[i]^Integer(v[i])

        if x==h:
            return vector(GF(5),v)

    raise ValueError("not in span")



#############################################################
# Matrix of tau on Cl(L_d)[5]
#############################################################

T=Matrix(GF(5),r5,r5)


for j in range(r5):

    I=basis[j].ideal()

    J=tau(I)

    h=Cl(J)

    T.set_column(j,coordinates(h,basis))



print("tau matrix:")
print(T)



#############################################################
# Sanity checks
#############################################################

print("tau^4 = 1 ?")
print(T^4 == identity_matrix(GF(5),r5))


#############################################################
# chi_d component
#############################################################

dim=(T-2*identity_matrix(GF(5),r5)).right_kernel().dimension()


print("--------------------------------")
print("dim Cl(L_d)(chi_d) =",dim)
print("--------------------------------")
