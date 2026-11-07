# nil31 short r=2/3: a fixed modular candidate cannot lift

2026-11-07. Task `E_multiple_hyperplane`. This note eliminates one specific
untried construction, and proves no component membership or global exhaustion.
The canonical full source and its dense single-orbit reduction are
`MODULAR_NIL31_SHORT_SUPPORT.md` §7. We retain its full t=d=1 tensor, including
its constant terms.

Put s=x, u=w+εx and qε=y−zu+ε(y²−u²)/2. On (y,z,u), define

    ιπε(dy∧dz∧du) = (2z dqε − qε dz)/3.

Its three brackets are

    {y,z} = −2z(z+εu)/3,
    {y,u} = [y+zu+ε(y²−u²)/2]/3,
    {z,u} = 2z(1+εy)/3.

These satisfy the full three-variable Jacobi identity identically for every ε.
Their modular field is (z+εu)∂y+(1+εy)∂u. If the complete four-variable
extension could retain −div V=y, this would give exactly the rank-two
semisimple modular-field candidate suggested in §7, after returning to x,w.
The required initial derivation is

    V0 = (z−y+2zu−y²/3)∂y + (z−2yz/3)∂z + (1+yu/3)∂u.

It does preserve π0. However, there is **no regular formal extension**
Vε=V0+εV1+… through quadratic x-independent Poisson derivations of this fixed
πε, even if the condition on div V is completely omitted.

Indeed let a be the coefficient of yu in V1_z. In the yz component of
L_(V0+εV1)(π0+επ1), the εy² coefficient is −a/3 and the εyzu coefficient is
−5a/3−2/3. They cannot both vanish. Equivalently, −5 times the first coefficient
plus the second is the nonzero constant −2/3. This is already a first-order
contradiction, and cannot be repaired by any higher term of Vε.

The standalone exact producer
`code/classification_nil31_short_semisimple_candidate.py` reconstructs the whole
48×30 first-order system over QQ and this two-row left certificate. Its matrix
has rank24 and its augmentation rank25. The report is
`reports/CLASSIFICATION_NIL31_SHORT_SEMISIMPLE_CANDIDATE_20261107.json`; the finite
contract is `jobs/nil31_short_two_thirds_semisimple_deformation.json`.
Measured runtime was 0.171 seconds. Preliminary exact degree-one and degree-two
polynomial ansätze were inconsistent as well; the first-order certificate is
the sufficient result, so no higher-jet search is warranted here.

Geometrically, the chosen qε smooths the restriction to the h=z=0 plane: the
original two hyperplanes belong to the same tangent ruling of the homogeneous
quadric, whereas generic ε breaks that incidence. A valid deformation must
extend the entire Poisson derivation; integrability of πε and semisimplicity
of the proposed modular field alone do not provide it.

This excludes only this fixed πε construction. It leaves open x-dependent
four-variable families, different translation directions, different three-
variable tensors, and all component membership questions for the original
source. The integrator independently reconstructed the two displayed coefficients.
For L_(V1)π0 in the yz position, the only terms contributing y² or yzu are
−4z(V1_z)/3 and −(y+zu)∂u(V1_z)/3. They contribute −a/3 and −5a/3.
For L_(V0)π1, the yzu contributions are 2/9, −4/9, −4/9, totaling −2/3;
the y² contribution is zero. This proves the two-row contradiction without
assuming the producer's full matrix rank or importing its code. The independent
review does not expand the scope to x-dependent or different Poisson families.
