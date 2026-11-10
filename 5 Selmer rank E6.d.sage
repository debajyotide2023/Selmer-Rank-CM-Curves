reset()

K.<a> = QuadraticField(-7)

print(K)

############################################################
# Elliptic curve
############################################################

d = 2

A = -35*d^2
B = -98*d^3

E = EllipticCurve(K, [0, 0, 0, A, B])

print(E)

############################################################
# 5-division polynomial
############################################################

R.<x> = PolynomialRing(K)

psi5 = E.division_polynomial(5)

print("psi5 factorization:")
print(psi5.factor())

############################################################
# choose factor (5-torsion x-coordinate factor)
############################################################

f = psi5.factor()[0][0]

print("\nChosen factor:")
print(f)

############################################################
# F = K(beta)
############################################################

F.<beta> = K.extension(f)

print("\n[F:K] =", F.relative_degree())
print("[F:Q] =", F.absolute_degree())

############################################################
# Verify beta is a 5-torsion x-coordinate
############################################################

print("\npsi5(beta) =")
print(beta.minpoly()(beta))
print(f(beta))

############################################################
# L = F(gamma)
############################################################

S.<Y> = PolynomialRing(F)

gamma2 = beta^3 + A*beta + B

L.<gamma> = F.extension(Y^2 - gamma2)

print(L)

print("\n[L:F] =", L.relative_degree())
print("[L:Q] =", L.absolute_degree())

############################################################
# Check relation
############################################################

print("\ngamma^2 =")
print(gamma^2)

print("\nExpected =")
print(gamma2)

print("\nRelation correct?")
print(gamma^2 == gamma2)

############################################################
# Point P
############################################################

P = E(L)([beta, gamma])

print("\nP =")
print(P)

print("\n5P =")
print(5*P)

############################################################
# CM element
############################################################

omega = (1+a)/2

print("omega =")
print(omega)

print("\nomega^2-omega+2 =")
print(omega^2-omega+2)

############################################################
# Silverman's formula
############################################################

D = 7*(1-omega)^4

beta0 = beta/d

beta0_1 = (beta0 - D/(beta0-omega^2+2))/omega^2

beta1 = d*beta0_1

print("\nbeta1 =")
print(beta1)

############################################################
# Verify beta1 is again a root of psi5
############################################################

print("\npsi5(beta1)=")
print(f(beta1))

############################################################
# theta from Silverman's y-coordinate formula
############################################################

theta = (1/omega^3)*(1 + D/(beta0-omega^2+2)^2)

print("\ntheta =")
print(theta)

print("\nCheck theta^2:")
print(theta^2 == (beta1^3 + A*beta1 + B)/(beta^3 + A*beta + B))

############################################################
# New point P1
############################################################

P1 = E(L)([beta1, theta*gamma])

print("\nP1 =")
print(P1)

print("\n5P1 =")
print(5*P1)

#############################################################
# Define CM automorphism tau(P) = omega*P
#############################################################

Labs = L.absolute_field('z')

phi = L.embeddings(Labs)[0]

a_abs = phi(a)
beta_abs = phi(beta)
gamma_abs = phi(gamma)

target_beta = phi(beta1)
target_gamma = phi(theta*gamma)

Aut = Labs.automorphisms()

tau = None

for s in Aut:
    # tau must lie in Gal(L_d/K), i.e., it fixes a = sqrt{-7}
    if s(a_abs) == a_abs:
        if s(beta_abs) == target_beta and s(gamma_abs) == target_gamma:
            tau = s
            break

if tau is None:
    raise RuntimeError("CM normalized tau not found")

print("tau found:")
print(tau)

print("\nCheck tau fixes sqrt{-7}:")
print(tau(a_abs) == a_abs)

print("tau(a) =")
print(tau(a_abs))

print("a =")
print(a_abs)

#############################################################
# Class group
#############################################################

Cl = Labs.class_group(proof=False)

print("Class group:")
print(Cl)

#############################################################
# Extract 5-torsion
#############################################################

five_candidates = []

for g in Cl.gens():
    n = g.order()
    if n % 5 == 0:
        h = g^(n//5)
        if h.order() == 5:
            five_candidates.append(h)

def in_span(h, basis):
    r = len(basis)
    for v in GF(5)^r:
        x = Cl(1)
        for i in range(r):
            x *= basis[i]^Integer(v[i])
        if x == h:
            return True
    return False

basis = []

for g in five_candidates:
    if not in_span(g, basis):
        basis.append(g)

r5 = len(basis)

print("5-rank =", r5)

print("Basis of Cl(L)[5]:")
for b in basis:
    print(b)

#############################################################
# Coordinates in the 5-torsion basis
#############################################################

def coordinates(h, basis):
    r = len(basis)
    for v in GF(5)^r:
        x = Cl(1)
        for i in range(r):
            x *= basis[i]^Integer(v[i])
        if x == h:
            return vector(GF(5), v)
    raise ValueError("Element not in span of basis")    

#############################################################
# Matrix of tau on Cl(L)[5]
#############################################################

T = Matrix(GF(5), r5, r5)

for j in range(r5):
    I = basis[j].ideal()
    J = tau(I)
    h = Cl(J)
    T.set_column(j, coordinates(h, basis))

print("tau matrix:")
print(T)

#############################################################
# Characteristic polynomial mod 5
#############################################################

R.<x> = PolynomialRing(GF(5))

charpoly = T.charpoly()

print("Characteristic polynomial:")
print(charpoly)

# Minimal polynomial of omega mod 5 is x^2 - x + 2 (since 5 is inert in K)
m = x^2 - x + 2

fac = charpoly.factor()
print(fac)

e = charpoly.valuation(m)

print("Multiplicity of omega's minimal polynomial =", e)
