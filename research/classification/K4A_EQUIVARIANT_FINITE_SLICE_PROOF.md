# K4a：C* weight で有限となる exact transverse slice の証明

2026-11-08、独立幾何 audit。
本稿は有限 producer 証明書の**論理**を証明する。新しい点、jet、minor 又は nilpotence 指数を要求しない。
具体的 polynomial graph／全三百五十式の結果は、別の immutable producer と独立 arithmetic receipt で確認する。
結果がまだ無い段階で二十三／二十四のどちらの component も断言しない。

## 1. raw 空間と実際の群

V は五変数 homogeneous quadratic bivector の百五十係数空間。
Z⊂V は元全三百五十 projective Jacobi equations `E∧[P,P]=0` の scheme。
Euler redundancy `P→P+E∧L` は L linear modulo C*E の二十四次元 additive group。
homogeneous weight零から、これは全 equation を正確に保つ algebraic translation である。
PGL5 と nonzero tensor scale も Z を保つ。
従って実際の smooth algebraic group

`H=A24⋊(PGL5×Gm)`

が作用する。GL5 の scalar center は homogeneous quadratic bivector へ自明に作用するので二重に数えない。
実 projective tensor が非零の近傍では Euler quotient 二十四と tensor projectivization 一は free、
raw dimension から二十五を引いて projective tensor の dimension を得る。

保存済み K4a point P は元全 Jacobi を満たす。
既存の独立 certificate は `rank D=101`、実 H-orbit tangent `dim G=48` を証明する。
従って kerD は四十九次元。
保存済み N は DN零、既存 all51 tangent columns の annihilator が N で8/87なので N∉G。
新しい百五十列 minor を作らず、`kerD=G⊕C N` が従う。

## 2. 同じ normalized torus で全 weight を読む

coordinate order(c,a,b,u,v)、weights(2,1,1,0,0)を用い、実作用を

`λ^-1 diag(λ²,λ,λ,1,1)_*`

とする。これは P を固定する stabilizer torus。
raw coefficient indexed by pair(i,j), degree-two monomial e の character は

`wi+wj−Σek wk−1`。

全値は[−5,2]、保存済み N はこの**同じ**作用で weight−1。
tensor rescaling を入れる前の weight と混同しない。
degree-four equation row indexed by quad I, monomial e の character は

`Σi∈I wi−Σek wk−2`。

従って各元 row は固有 equation。これは full polynomial equation の covariance であり、
tangent matrix の support のみから非線形 equivariance を仮定したものではない。

## 3. actual linear slice と selected implicit equations

既存の nonzero百一Jacobi minor の coefficient columns を U とする。
各 column は raw coordinate vector、従って torus eigenvector。
D|U 可逆なので V=kerD⊕U。
slice は実際の affine linear subspace

`S=P+U+C N`、dim S=102。

N∉G から S の tangent は G の直和 complement。
P が固定され、U/N が eigen なので S は torus-stable。
既存 minor の百一 equation rows を F1,…,F101 とする。
slice coordinates を x∈U,t∈C N とすると ∂F/∂x(P) が可逆。
formal implicit function theorem は唯一の graph

`x_i=f_i(t), f_i(0)=0, f_i'(0)=0`

を与える。最後の零は DN=0 から従う。
この段は selected equations の formal **全解**であって有限 truncation の existence 仮定ではない。

## 4. weight が graph の全 tail を消す

selected eigen equations と slice は equivariant、implicit graph は唯一なので equivariant。
weight(x_i)=χi、weight(t)=−1 とすると

`f_i(λ^-1t)=λ^χi f_i(t)`。

formal power series の各 monomial t^j を比較すると、非零項は j=−χi の一つだけ。
χi∈[−5,2]、constant/linear terms 零なので可能なのは j=2,3,4,5。
従って **全 formal graph は degree五以下の literal polynomial graph**。
大きい次数の coefficient を未計算のまま零とみなしたのではない。

有限 producer は weight-blockごと j2,…,5 の係数を一括決め、
actual family `P+tN+Σj=2..5 t^jTj` を出せばよい。
各 Tj は U 内の weight−j 部分のみ。
元 equations は quadratic なので family 上の全 row は degree十以下の polynomial。
selected百一と full三百五十の全 degree0,…,10 coefficient を正確に代入確認する。
これで selected graph の存在と全式の残余が有限にレビュー出来る。
局所化した formal jet を将来伸ばす義務は生じない。

## 5. 全 row の零／非零による reduced support

selected graph 上の full equations を gj(t) とする。
full completed slice ring は正確に `C[[t]]/(gj(t))`。

* 全 gj が零なら completed slice は smooth formal curve。
  literal polynomial graph 自体が actual algebraic curve を与える。
* 一つでも非零なら一変数 formal DVR で radical ideal は(t)。
  completed reduced slice は一点。
  exact exponent は reduced component の義務に不要。

N が tangent normal であることだけから第二の場合を排除しない。
逆に finite polynomial の実代入が零であることを確認した場合だけ第一の場合を主張する。

## 6. étale action product：closed orbit の仮定は不要

H の identity 近傍で orbit derivative が G の basis となる四十八 local parameters を取る。
例えばその四十八 tangent directions を持つ smooth local subvariety T⊂H を選べばよい。
T は subgroup である必要はない。

`T×S→V, (h,s)↦h*s`

の differential は G⊕TS=V のため同型。従ってこの map は(P,id)でétale、
completed local ambient rings は同型となる。
H は Z を実際に scheme として保つので preimage は T×(Z∩S)。
ゆえに completed Z は smooth四十八変数と completed slice の product。
Luna slice 又は orbit が閉じているという仮定を使わない。

非零 gj の場合は raw reduced local dimension四十八、
Pの orbit closure は四十八次元なので full irreducible component。
実 tensor projectivization 後は二十三次元となる。
全 gj 零の場合は raw smooth local dimension四十九、projective二十四次元。
唯一の local component がその actual polynomial curve と group action で現れる。

orbit は自分の closure 内で open なので、normal derivative N を持つ actual curve が
同じ二十三次元 orbit closure 内にあるとは言えない。
全零の場合、stabilizer torus は t を nonzero scalar 倍するため全 t非零は同じ実 projective orbit。
これは E親認識に利用可能だが、どの既 component か又は新 component かは別の distinctness 義務。

## 7. 独立 certificate の必要チェック

独立側は既存 D101/G48/fullkernel 証明書を再利用し、新しい profile を増やさない。
P/N/U/selected rows の既存 indices、normal functional8/87、同じ normalized torus、
graph の support/weight、selected implicit equations、full三百五十 row の全 polynomial coefficients を確認する。
producer の一次行列又は SymPy code を import して結果を繰り返すだけにはしない。
数学的二分岐のうち、実際の full coefficient receipt が証明した側だけを採用する。
