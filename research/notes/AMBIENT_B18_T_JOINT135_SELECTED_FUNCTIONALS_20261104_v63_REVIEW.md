# Selected joint135 third-order functionals: independent review

The primary computes new restricted third-order implicit functionals, not
a general105-source or30-gauge third jet. This review uses the independent
sparse trace algebra, physical tensor push and flag/frame implementation.
It postverifies the given data without solving for any source or gauge
coefficient.

The original source coefficient derivatives are rebuilt by independent
sparse AD. Each of the135 physical/gauge input directions is then evaluated
in its source-dependent quadratic trace algebra. All30 canonical constraint
and34 G8/nu target derivative rows are reproduced. In particular the monic
relations vary with the physical source, and the nonzero odd gauge center
is retained. The scalar block Jacobian is the actual
[[J105,0],[J_gauge_physical,J30]].

A separate25-monomial source-DAG replay computes the constraint residuals
with all unknown cubic coefficients zero. The tuple permutation between
the earlier joint inventory and this primary is explicit. All105 source
equations and30 canonical scalar constraints agree on the known base20.
Their five cubic residuals, and every entry of the nonlinear known G8/nu
response, are verified. The old physical105 pure-r³ coefficients supply
an independent source crosscheck; no previously complete gauge r³ cache
is assumed.

The entire native source-origin vector is checked to vanish for every one
of619 coherent combinations before omitting the N3 term. All34 primitive
G/nu derivatives and all requested native covectors are reproduced using
the actual N0 and scalar monic reductions. Each of152 literal adjoints is
then checked against all135 Jacobian components. Their scalar products
with the actual residuals reproduce all619 by five implicit contributions.

These last values are only -rho*f. The direct nonlinear lower response q
is independently checked in the earlier phase and must be included exactly
once in the later source/native assembly. Thus the selected-method formula
is q-rho*f; this review does not turn the implicit part alone into a complete
nonlinear target value. No full619 quartic value, actualHF4, or exact totalT
length is asserted.
