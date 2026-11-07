# rank零plane・m二・半単純normal型のCartan枝は全MPS40へ入る

2026-11-07、主義務G。producerによる一次scope確認と全source幾何の橋。
独立review前はglobal採録を待つ。

**命題。** rank零planeに沿うm_D=2 first symbolのtraceless co-Higgs部分が
一つのglobal vector field Vとconstant semisimple sl2 matrixの積であり、
trace field TとVがsimultaneously diagonalizableとする。
このfirst symbolを持つ全higher Poisson tensorは、既知MPS40 componentの和に属する。
generic eigenvalue、nonresonant period、higher-weight消去を追加で仮定しない。

## 1. whole weight零sourceと実SNC point

normal GL2でconstant semisimple matrixをdiagonalに直す。
traceless部分の非零scaleをVへ吸収すると、全weight零sourceは

`Π0=zX∧∂z+wY∧∂w+K(z,w)∂z∧∂w`、
`X=T+V`, `Y=T−V`

である。
全weight零Jacobiは `[X,Y]=0` と `K=bzw` を与える。
実際base関数f,z,wのJacobiは
`zw[X,Y]f+(az²−cw²)(X−Y)f`
（K=az²+bzw+cw²）であり、X−Y非零からa=c=0、[X,Y]=0。
projective Euler wedgeも同じlinear-field条件を与え、trace-free liftのcommutatorが
scalar Eulerとなる場合はtraceからそのscalarも零である。

simultaneous diagonalizationによりbase座標x0,x1,x2を取り
`X=Σα_i x_i∂i`, `Y=Σβ_i x_i∂i` とする。
m_D=2は `det(Euler_base,X,Y)` が恒等零でないことと同値。
従って

`κ=det[[1,1,1],[α0,α1,α2],[β0,β1,β2]]≠0`。

actual normal-weight零Pfaffianは符号を除いて

`D0=κ*x0*x1*x2*z*w`。

これは五つの独立hyperplaneのsimple normal crossings divisorである。
Π0はtoric logarithmic tangent bundle上のconstant skew matrixで、Pfaffianκは非零なので
全域でlog symplecticである。
b又は一部biresidueが零でも、κ非零という全条件を保つ限りこの結論は変わらない。

## 2. 一次資料の正確な全SNC scope

M. Matviichuk–B. Pym–T. Schedler,
*A local Torelli theorem for log symplectic manifolds*, arXiv:2010.08692v2,
<https://arxiv.org/html/2010.08692v2>、及び
<https://arxiv.org/pdf/2010.08692v2> を直接読解した。
Theorem1.4（印刷7頁）はP4のPoisson moduliのうちSNC log-symplectic locusと交わる
irreducible componentを正確に四十個と分類する。
§5.3（印刷39–41頁）とExample6.24（印刷55–56頁）がその全scopeを支える。

特にExample6.24はprojective toric boundaryでresonant local-system obstructionが消えるため、
全toric log-symplectic periodでformal neighborhoodをsmoothing diagramsにより記述する。
n=2の四十familiesのglobal distinctnessを直接確認しTheorem1.4へ結ぶ。
Theorem1.2の一般のnonresonant openだけを全SNC periodへ拡張したものではない。
同論文が未解決としているのは、全log-symplectic componentがSNC locusに達するかどうかであり、
本命題では次節でその実際の到達を別に証明する。

P4のrigidityとconnected PGL5 actionにより、pair deformation又はPGL quotientとしてのmoduliの
irreducible componentは、固定P4上のPoisson parameter spaceのcomponentと対応する。
したがって既存catalogueのMPS40へ適用できる。
外部定理から全Poisson componentが四十個だとは主張しない。

## 3. 全higher sourceからのactual C* degeneration

conormal abelian性によりnegative normal-weight項は無く、全quadratic bivectorは

`Π=Π0+Π1+Π2`

（全weight一30、全weight二9）と書ける。
actual PGL5 one-parameter subgroup
`g_t=diag(1,1,1,t^−1,t^−1)` のpushforwardは、t≠0で

`Π_t=Π0+tΠ1+t²Π2`

となる。
これはformal jetではなく全quadratic tensorのactual algebraic familyであり、
元点t=1とSNC点t=0を含む。
Jacobiはtについてpolynomialなので全tで成立する。

connected algebraic groupは有限型varietyの各irreducible componentを保存する。
例えばC*×Cのimage closureはirreducibleでCを含むため、maximalityからC自身に入る。
従って元Πを含む任意のglobal component Cは、全t≠0のΠ_tとそのlimitΠ0を含む。
Π0は実SNC log-symplectic pointなので、CはTheorem1.4の四十componentの一つである。
元Πもその和へ入る。命題が従う。

ここでは単に「limitΠ0が既知componentのboundaryに入る」ことから結論したのではない。
limitそのものがSNCであり、**それを通る全global component**が外部分類の射程に入ることを使う。

## 4. 残すwhole first-jet strata

constant normal sl2 matrixがnilpotentの場合、又はT,Vがcommuteしてもsimultaneously diagonalizableでない場合、
D0はnonreducedとなり得る。
そのΠ0がSNC tensorからさらに退化して得られるだけでは、本命題のcomponent argumentは使えない。
これらの全higher sourceには別の全30/9必要operator又は実component認識の証明が必要である。
またrankM一のO(1) co-Higgs部分も、一つのglobal Vのsupersourceとして保持する。
plane m_D≥3はnegative normal weightを持つ可能性があり、本C* familyを流用しない。
