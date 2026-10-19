# Enumeration algorithm

Let $M_{n,k}$ and $I_{n,k}$ count all and irreducible open singular meanders with $n$ transverse crossings and $k$ tangencies. The requested region is $n+2k\le W$, inclusively, for a nonnegative cutoff $W$. Write $M(x,t)=\sum M_{n,k}x^nt^k$ and similarly for $I$. The empty object is excluded: $M_{0,0}=0$; the two one-contact types contribute $x+t$.

## Frontier enumeration

[The transfer](Source%20files/irreducible_meander_enumeration/enumerate_irreducible_meanders.hpp) counts connected closed singular meanders $C_{n,k}$ with positive even $n$. Here the road is a simple closed curve and the river is an arc with distinguished endpoints. Equivalence preserves the disk orientation and both river endpoints, as in [Singular Meanders](https://arxiv.org/abs/2512.23091).

Sweep contact points along the river. A state records the noncrossing pairing of exposed ends ordered along the sweep cut, the split separating the ends below and above the river, and $k$. The pairing is encoded as a balanced bracket word. At contact layer $v=n+k$, the transverse count is $n=v-k$ and the weight is $v+k$.

Start with the empty frontier of multiplicity one. At each crossing or tangency, insert a pair of ends, continue a strand, or join adjacent ends and reconnect their former partners. The two road half-edges at a contact lie on opposite sides of the river for a crossing and on the same side for a tangency. Joining partners would complete a component and is rejected unless this closes the whole frontier; record that closure without extending it further. Distinct local configurations contribute separately even when their successor states coincide.

Equal states have identical continuations, so their multiplicities are added. A state and its reflection use the same canonical key, also summing multiplicities; final counts are not divided by reflection. Discard states whose remaining weight or crossing budget cannot close the exposed ends.

## Recovery of open counts

Exact generating functions supply rows $M_0,\ldots,M_4$, where $M_n(t)=\sum_kM_{n,k}t^k$. In particular,

$$M_0(t)=\frac{t}{1-t},\qquad M_1(t)=\frac1{(1-2t)\sqrt{1-4t}},\qquad M_2=M_1'/4,\quad M_4=M_3'/8.$$

Primes denote differentiation in $t$. The complete fixed-row expressions are in `all_meander_rows` in [meander_series.hpp](Source%20files/irreducible_meander_enumeration/meander_series.hpp). The transfer therefore records only completions with $n\ge6$; partial states with fewer crossings remain if they can reach such a completion. For $n\ge5$, recover the remaining open counts by the exact rank identities

$$M_{n,k}=\begin{cases}
\dfrac{n+1}{n+k+1}C_{n+1,k},&n\text{ odd},\\[4pt]
\dfrac{k+1}{2(n+k+1)}C_{n,k+1},&n\text{ even}.
\end{cases}$$

For $W\ge14$, enumerate through transverse cutoff $R=W-8-(W\bmod2)$ and closed weight $W_c=W+2-(W\bmod2)$.

## Prime factorization

Let $S$ count iterated snakes with $x+t$ removed. Its cubic equation determines coefficients recursively; maintaining $S^2$ avoids repeated cubic convolutions. Subscripts `odd` and `even` below denote parity in $x$.

For $P=x+O((x,t)^2)$ and $Q=t+O((x,t)^2)$, let $F_{P,Q}$ be the unique series satisfying $F_{P,Q}(P,Q)=x+t$. [Canonical prime factorization](https://arxiv.org/abs/2112.10289) gives

$$I=F_{x+S_{\rm odd},\,t+S_{\rm even}}-F_{M_{\rm odd},\,M_{\rm even}}.$$

Compute each inverse by eliminating coefficients in increasing ordinary degree $n+k$, where substitution is unit triangular. Weighted degree alone is insufficient: $Q$ can contain $x^2$, of the same weight as $t$. Here $P$ is divisible by $x$ and has weight at least one, while $Q$ has weight at least two. Substitution therefore cannot lower an outer monomial's transverse degree or weight: terms with $n>R$ or $n+2k>W$ cannot affect the retained coefficients.

The [forward calculator](Source%20files/meander_calculation/meander_calculator.hpp) first computes iterated irreducibles as the fixed point $J=I(x+J_{\rm odd},t+J_{\rm even})$. It then solves

$$A=S(x+B_{\rm odd},t+B_{\rm even}),\qquad B=J(x+A_{\rm odd},t+A_{\rm even})$$

by alternating substitutions, starting with $A=S$ and $B=J(x+A_{\rm odd},t+A_{\rm even})$, and returns $M=x+t+A+B$ when the truncated series stabilize. Completeness is limited by the supplied irreducible data.

## Boundary coefficients

Use $I_{0,k}=0$ and $I_{n,k}=0$ for $k<3$. Exact rows $I_1,\ldots,I_4$ are supplied by `irreducible_meander_rows`. With Euler's totient $\varphi$, for $m\ge1$,

$$I_{2m,3}=\varphi(m+4)-2,\qquad I_{2m-1,3}=0,\qquad I_{2m-1,4}=m\bigl(\varphi(m+4)-2\bigr).$$

For the remaining odd rows, cyclic rerooting gives

$$I_{n,k}=\frac{2(n+1)}k I_{n+1,k-1}\qquad(n\ge5\text{ odd},\ k\ge5).$$

Every remaining even coefficient has $k\ge4$, hence $n\le R$. Each odd coefficient above uses an even source with weight one less and transverse degree at most $R$. These cases cover the full requested triangle. Cutoffs below 14 require formulas only. Formula and cyclic values are checked wherever inversion independently supplies the same coefficient.

## Storage

[Frontier storage](Source%20files/irreducible_meander_enumeration/frontier_storage.hpp) uses full-key comparison in sharded hash tables. Transfer multiplicities use checked 128-bit integers; signed arbitrary-precision integers handle inversion. Divisions in rank identities are checked for integrality. The forward calculator uses the chosen `number_type`, which must provide sufficient precision.

If destination tables reach the memory budget, write their accumulated counts and subsequent transitions to hash buckets. Reduce these records separately, splitting oversized buckets without extending source states again or recounting completions. Completed layers replace their source files. Temporary files belong to a unique run directory and are removed when the call ends.
