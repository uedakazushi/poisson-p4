# Three exact homogeneous degree-eight seed identities

For the same literal homogeneous24-by56 raw matrix M^h used in the degree7
seed certificate, three of the112 reconstructed candidate columns have degree8:
zero-based candidate indices81,106,109. All three belong to the original QQ
module. The saved56-by3 matrix of linear homogeneous cofactors satisfies
M^h L=S_8 exactly. The native witness is a224-by3 rational coefficient matrix,
with row index4*j+k for original column j and ordered cofactor monomials
(w,P,e,a). Its complete size is1.9KB.

Expanding into all165 degree8 monomials and24 components gives complete
coefficient matrices C of shape3960-by224 and B of shape3960-by3. A field pivot
search chooses an invertible224-by224 square. Exact rational Dixon solving and
full QQ multiplication verify ALL3960-by3 coefficient identities, in0.044seconds.
The final proof is the literal QQ identity, independently of the modular search.

Together with the separate44 degree7 identities,47 reconstructed seeds are now
certified redundant QQ generators. The remaining65 candidate columns are not
proved to belong. These identities do not prove a QQ GB, localized surjectivity,
original product rank41, full cubic rank57, or an ambient statement. The distinct
augmented QQ GB attempt hit its60-second cap before producing a basis; its
failure does not affect these exact identities. This packet awaits independent
review.
