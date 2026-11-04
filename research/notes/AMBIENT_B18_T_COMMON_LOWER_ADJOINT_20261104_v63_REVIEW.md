# Independent common lower-adjoint/source-value review

This review checks the121 shared lambda functionals and their actual cached
values, the common106-entry covector D, all119 internal kernel checks, and
its619 source contractions. The619 lower values are not evaluated: nine zz
pairs and additional typed-root second coefficients are outside the previous
B11/762-root cache. No rawP4, X3 or higher source value is computed.

The checker reconstructs the84 unitstar and37 selectedL source functionals
from original768 labels and its own accepted Gamma4. It evaluates their
typed old B roots forwards in B11, retaining absent atom second coefficients
as symbolic variables. It independently replays the common105/copiedCB57
source through degree3 for the rawP3 terms, keeping regular quotient second
and third coefficients symbolic. All unavailable sensitivities in these
shared lambda values cancel before acceptance. Native effective S2 is used
only for the separately verified zero-origin shared source columns.

For the lower source z, the exact equations are

    b=B*C2, B=U0^(-1)*kappa2*m0,
    l=L*(pi-Ystar1*B)*C2,
    a=A*(I-Xstar1*B-XL1*L*(pi-Ystar1*B))*C2.

Here A is the fixed27 pivot left inverse and L is the fixed177-dimensional
lift through the37 linear generators times variables. Pulling
lambda_a*a+lambda_b*b+lambda_l*l back gives

    D=lambda_a*A+(lambda_l-lambda_a*A*XL1)*L*pi
      +(lambda_b-lambda_a*A*Xstar1
        -(lambda_l-lambda_a*A*XL1)*L*Ystar1)*B.

The independent checker computes A by inverting its actual selected27-row
unit minor, instead of replaying the primary's individual projection columns.
After the a correction, it verifies all119 linear-generator syzygies are
annihilated. It computes the177-space dual by descending triangular
coefficient substitution with the three free dual coordinates zero. It then
applies the output-Euler contraction forwards. All intermediate lambda,
projection and dual arrays and the final106 entries of D agree exactly with
the primary.

Every619 contraction D(Xc2-Xp1*a1) is independently folded to actual source
labels and checked, including the negative first-correction term. The nine
additional zz pairs are recorded as missing inputs, not assigned zero. The
native part continues to use complete zero-origin source combinations;
individual native second coefficients are not asserted equal to the
effective array. The result is a valid common covector and a concrete
source-availability frontier, not complete619 quartic values or HF4.
