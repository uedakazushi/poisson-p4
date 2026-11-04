# Admissible adapter changes are source equivalences

This packet proves an exact source-equivalence lemma and records a necessary
limitation on finite-jet admissibility. It does not reconstruct N2, prove that
an arbitrary free second coefficient is a gauge change, or decide HF4.

## 1. The ring and the exact admissibility equations

Work first in the ORIGINAL regular coefficient ring C9. The accepted actual
proper-core membership and source-interface packets give a regular sequence
G=(G1,...,G8), a unit12-by12 adapter N and the same physical column F12 with

    N F12=d*(G,0,0,0,0)^t.

The native syzygy matrix is N^t S(G), where S(G) has the28 Koszul columns on
G in its first eight rows and four free unit columns in its last four rows.
The old physical Hadamard columns B do not depend on N.

Let N'=E N be another exact adapter for the same physical F12 and the same
divisor d, with E near the identity. Write E in8+4 blocks as

    E=[[A,B],[C,D]].

Exact admissibility means G'=A G and C G=0. The latter is an all-order
equation in C9. Regularity of G supplies a regular matrix H with

    C^t=K(G)*H,

where K(G)[:,i<j]=Gi*ej-Gj*ei. No regularity of G after ell or in an Artin
algebra is used.

## 2. The literal source-unit matrix

Koszul naturality with this sign convention gives

    A^t K(A G)=K(G)*wedge2(A^t).

Consequently

    E^t S(A G)=S(G)*T,
    T=[[wedge2(A^t),H],[B^t K(A G),D^t]].

The checker verifies every224 naturality entry and every384 entry of this
identity as universal polynomials in all entries of A,B,D,H and G. The
identity does not use special numerical matrix values or a source point.

At the C9 closed point, G=0 and A=D=I. Thus T is upper triangular with
identity diagonal, regardless of H at that point. It is a unit32-by32 matrix.
Multiplying the identity by N^t gives

    N'^t S(G')=N^t S(G)*T.

Restriction of scalars in the actual four-element monic basis gives a
unit128-by128 source matrix. The full comparison X changes by
diag(I640,Res4(T)); its old physical640 columns remain literally the same.
Subsequent ell or Artin base change preserves this matrix identity and its
inverse. This is an exact coherent SOURCE change, not a common target-row
change involving the old B matrix.

Thus admissible adapters give the same complete presentation and intrinsic
initial Hilbert function. A particular fixed619 coordinate vector need not
remain identical without transporting its source coefficients and lower
Schur corrections. Partial-family numerical invariance is not inferred from
the complete source-unit statement alone.

## 3. What finite-jet admissibility actually requires

The exact parameterization is

    E=[[A,B],[(K(G)H)^t,D]],  G'=A G,

subject to E being a unit and any prescribed order of agreement. The
admissible jets are the image of these all-order parameter jets. They need
not equal all free adapter jets satisfying only a truncated annihilator
equation. Such a replacement would require a separate lifting proof.

For example, let R=k[[x,s]], C=R[t]/(t^2-s), with base ideal m=(x,s), and
G=(x,t), a regular sequence in C. Modulo m^3 the row (0,x^2*t) annihilates G,
because x^2*t^2=x^2*s is zero there. But it is not the reduction of an exact
syzygy. Every exact syzygy is (-t*h,x*h); writing h=a+b*t, the functional

    coefficient(first even,x*s)+coefficient(second odd,x^2)

is -b_x+b_x=0. Its value on (0,x^2*t) is1. The checker verifies this on all12
monomial h basis vectors of C/m^3. This is a formal algebra countercondition,
not a new physical/source point.

In particular, neither a free second-adapter perturbation nor the choice
N2=0 is declared admissible. The current actual reconstruction should use
the existing B11 physical/gauge/coordinate values and the full rational
telescopes, checking NF12=d(G,0) in the needed channels. The lemma guarantees
coherent invariance once exact admissibility is established; it does not
supply missing coefficient values or waive that source check.
