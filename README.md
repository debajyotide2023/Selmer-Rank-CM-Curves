# SageMath computations for two CM quadratic twist families

This repository accompanies the paper's algorithm for computing $5$-Selmer dimensions in the quadratic twist families of $X_{\mathrm{ns}}^+(11)$ and $X_0(49)$. The two scripts compute the class-group quantities used in the ordinary and supersingular cases of the algorithm.

| Family | CM field | Reduction at $5$ | Quantity printed |
| --- | --- | --- | --- |
| $E_{9,D}$, twists of $X_{\mathrm{ns}}^+(11)$ | $\mathbf Q(\sqrt{-11})$ | Good ordinary | Dimension of the $2$-eigenspace of $\tau$ on $\mathrm{Cl}(L_d)[5]$ |
| $E_{6,D}$, twists of $X_0(49)$ | $\mathbf Q(\sqrt{-7})$ | Good supersingular | Multiplicity $e$ of $x^2-x+2$ in the characteristic polynomial of $\tau$ on $\mathrm{Cl}(L_d)[5]$ |

Both quantities are intermediate inputs to the paper's $5$-Selmer computation.

## Run the examples

Open a terminal in the repository folder and run the scripts with SageMath.

```bash
sage "5 selmer rank E9.d.sage"
sage "5 Selmer rank E6.d.sage"
```

- [E9 ordinary SageMath script](./5%20selmer%20rank%20E9.d.sage)
- [E6 supersingular SageMath script](./5%20Selmer%20rank%20E6.d.sage)

## Ordinary case (E9)

The base model is

```math
E_9\colon y^2=x^3-9504x+365904.
```

The example sets `d = -14`, corresponding to fundamental discriminant $D=-56$. Over $K=\mathbf Q(\sqrt{-11})$, the script adjoins a root $\beta$ of the quadratic $X^2+b_1X+b_2$ specified in the code. It then adjoins $\gamma$ satisfying

```math
\gamma^2=d^3(\beta^3-9504\beta+365904).
```

The point on the twisted model is $P=(d\beta,\gamma)$. The script constructs $L_d=K(\beta,\gamma)$ from this point and selects the automorphism $\tau$ fixing $K$ that acts as $\tau(P)=2P$. It does this by requiring $\tau(\beta)=-b_1-\beta$ and $\tau(\gamma)=c\gamma$, using the prescribed formulas for the coordinates of $2P$. For this example the field degrees are $[F:K]=2$, $[L_d:F]=2$, and $[L_d:\mathbf Q]=8$.

The script builds the matrix $T$ of the forward ideal-class action $[I]\mapsto[\tau(I)]$ on $\mathrm{Cl}(L_d)[5]$. The line `dim Cl(L_d)(chi_d)` reports

```math
\dim_{\mathbb{F}_5}\ker(T-2I),
```

the dimension of the component for the character normalized by $\chi_d(\tau)=2$. The two values under `Check K embedding` should be `-11` and `0`. The gamma relation, the check that $\tau$ fixes $\sqrt{-11}$, and the matrix check $T^4=I$ should each print `True`.

## Supersingular case (E6)

The base model is

```math
E_6\colon y^2=x^3-35x-98.
```

The E6 script sets `d = 2` for the twist with fundamental discriminant $D=8$. Set `d = 1` to run the base model. Over $K=\mathbf Q(\sqrt{-7})$, the script constructs the $5$-division field by adjoining the coordinates of a nonzero point $P$ of order $5$. For $\omega=(1+\sqrt{-7})/2$, it selects the automorphism $\tau$ fixing $K$ that sends $P$ to $[\omega]P$.

For this example, the expected field degrees are $[F:K]=12$, $[L:F]=2$, and $[L:\mathbf Q]=48$. The two division-polynomial evaluations should print `0`, and both `5P` and `5P1` should print the point at infinity, `(0 : 1 : 0)`.

Here `Multiplicity e` is the exponent of the irreducible polynomial $m(x)=x^2-x+2$ in the characteristic polynomial of $T$ over $\mathbf F_5$. The action of $T$ has order dividing $24$, which is coprime to $5$. The corresponding component is $\ker(T^2-T+2I)$ and has $\mathbf F_5$-dimension $2e$.

## Twist parameter conventions

The paper indexes twists by fundamental discriminants $D$ with $5\nmid D$. In the scripts, the squarefree twist parameter is $d=D$ when $D\equiv1\pmod4$ and $d=D/4$ when $4\mid D$. The short model of a twist of $y^2=x^3+Ax+B$ is $y^2=x^3+Ad^2x+Bd^3$.

For the ordinary script, `d` is already an input near the top of the file. For the supersingular script, a change of `d` affects both the coefficients and the CM formula. For the $-98$ model, set $r=\omega^2-2$ and $C_\omega=7(1-\omega)^4$. The formulas on the twisted model are

```math
q=\beta-dr,\qquad
\beta_1=\frac{\beta-d^2C_\omega/q}{\omega^2},\qquad
\theta=\frac{1+d^2C_\omega/q^2}{\omega^3}.
```

The supersingular script calls the CM constant $C_\omega$ `D`. This is separate from the paper's fundamental discriminant $D$. For `d = 2`, the twisted curve is $y^2=x^3-140x-784$.

## Reproducibility and citation

The ordinary script calls `Labs.class_group()` and uses SageMath's current proof setting. The supersingular script explicitly uses `proof=False`, so its class-group result is conditional on GRH. To request unconditional class-group computations in both scripts, use `Labs.class_group(proof=True)`. This can increase the runtime substantially.

When reporting a result, record the family, the twist discriminant, the SageMath version, the class-group proof setting, and the repository release or commit. Cite the accompanying paper for the Selmer algorithm and the specific code version for the computation.
