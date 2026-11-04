# 低六不変式のnon-HSOP obstruction：独立レビュー

2026-11-04。**指定degree typeのHSOP不存在、実low6写像の非finite性、
商空間の正次元origin fiberを受理する。** NEWの定義式次数を決める結果ではない。
新しい点・不変式評価・source26・CAS・Sym冪DP再構築は行わない。

## 完全characterと整数convolution

対象はC上の `R=C[W16]^Sp4`。保存済み完全characterとその独立レビューを使用する。
後者は全16 weightsのWeyl DPを26次まで独立再構築し、Inv26全体の次元45を認証している。
構成した45本の独立性だけを上界と取り違えていない。
QQ上の完全character次元はchar0の係数拡大後も同じである。

今回の独立checkerはC2の四positive roots
`(2,0),(0,2),(1,1),(1,-1)` の全16 subsetからdenominatorの八符号付きshiftを作り、
全27 character recordのweight/signと整数和を照合する。重いDP自体は再実行しない。
必要な完全次元は

```
H2=H6=H10=0, H14=3, H18=9, H22=21, H26=45。
```

degree vector `(4,8,8,8,8,12)` の全64 subsetを直接数えて分母多項式を構成し、
逐次多項式積を使ったauthor実装とは別に26次係数を求める。

```
[t26]H_R(t)(1-t4)(1-t8)^4(1-t12)
=H26-H22-4H18+3H14+7H10-2H6-8H2
=45-21-36+9=-3。
```

未計算の高次数はこの係数に寄与せず、Hilbert級数全体の閉形式も不要である。

## CMからHSOP不存在へ

Sp4はchar0でlinearly reductive、C[W16]はregular polynomial ring。
Hochster–Roberts定理をこの標準的な状況に適用し、RはCohen–Macaulay。
Rは有限生成connected positive-graded domain、既受理の商関数体からdimR=6である。

このdegree typeのhomogeneous SOPが存在すれば、CM性により長さ6のregular sequenceとなる。
各元を順に割る短完全列から、そのquotientのHilbert級数は上記積になる。
Artinian quotientのHilbert係数は非負なので−3は矛盾。
したがって**いかなる係数選択でも、このdegree typeのHSOPは存在しない**。
これはCM性があるから言えることで、任意の非CM環へ同じ推論を移さない。

引用する標準定理はM. Hochster–J. L. Roberts,
*Rings of invariants of reductive groups acting on regular rings are Cohen-Macaulay*,
Adv. Math.13(1974),115–175, DOI10.1016/0001-8708(74)90067-X。
本監査はその証明や新規性を主張せず、仮定の適用を確認する。

## 実際の六つ、finite性、zero fiber

特定の `F4,p1,p2,p3,p4,N12` の独立性は別の受理済みglobal恒等式から従う。
`F4²=400p0+1600p4`、`N12=2p0*m` を微分し、
旧普遍rank6からの変数変換のdeterminantは `400p0/F4`。
同じ保存点でF4,p0は非零なので、この六つは代数的独立である。
六個あることやInv8の線形独立性だけから推定しない。

A=C[low6]は六変数weighted polynomial ring。FracR/FracAは同じ超越次数6なので有限、
したがってSpecR→SpecAはgenerically finite。
もしRがA上module-finiteなら六つがHSOPとなり矛盾するため、全域ではfiniteではない。

さらに `R/(low6)` が有限次元なら、その有限homogeneous basisをRへ持ち上げる。
任意のhomogeneous元をquotient basisで減じた残りはΣf_i*s_iであり、
deg(s_i)<deg(元)なので次数帰納により有限個の持ち上げがA-moduleとしてRを生成する。
**RのA上有限性を先に仮定するNakayamaではない。**
逆方向はmodule-finite商を取るだけ。
従ってこのquotientは無限次元。有限生成C-algebraのdimension0なら有限次元なので、
SpecR上のorigin fiberは正次元となる。

このfiberは全正次数不変式の共通零点である通常のnullconeとは異なる。
C上の非vertex quotient pointは存在するが、今回構成・評価していない。
必要ならReynoldsのR線型射影からaffine quotientの全射性も確認できる。
fiberの正確な次元や成分数は未決。

## Ordered chartと次数推論の限界

low6=0ならF4と全pencil係数が零。ordered chartでは
`F4=20(x²+y²),p4=x²y²` なので、**幾何的点として**x=y=0であり、xyDelta≠0の開から外れる。
これはideal `(x²+y²,x²y²)` と `(x,y)` のschemeとしての同一性を述べていない。
generic残余曲線で除外したboundaryを、originで無視できないことが分かる。

generic field degreeとNEWのgeneric multisection degreeだけを使うHilbert multiplicity計算は、
この否定されたfinite-HSOP仮定を補わない限り成立しない。
degree30の定義式の存在を証明したとも否定したとも言わない。
新しいHSOP、nonnull point、source26制限、NEW定義式はこのreviewの対象外。

再現:

```
python code/q_firstjet_low6_non_hsop_review_20261104_v63.py \
  --output /tmp/low6-non-hsop-independent-review
```
