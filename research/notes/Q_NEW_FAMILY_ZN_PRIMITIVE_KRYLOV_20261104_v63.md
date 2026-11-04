# ZNは10次shape拡大のprimitive element：同じ旧点による有限自由certificate

2026-11-04。**ZN=N₁₂/F₄³が既知sourceの10次shape体拡大を生成することを証明した。** 保存済みの同じβ点index0、同じ素数32003を一度だけ用い、有限自由10次algebraの係数行列式が非零になることを確認した。新しいβ点、関係式fit、巨大なN norm、source26、追加Inv₁₄展開は使っていない。

## 既に確定したfieldと線形s回収

K=QQ(R,Z₀,Z₁,Z₂)、L=QQ(R,Z₀,r,s)とする。raw determinant identity、full content最大性、二つのs写像の次数3,5はそれぞれ独立受理済みであり、保存したGはK[r]で既約な10次多項式。さらに

\[
 L=K(r),\qquad [L:K]=10.
\]

既存の正確なBrown PRS prefixのlinear polynomialを

\[
 C_1=A(r)s+B(r)
\]

と書く。実generic incidenceではA≠0で、s=−B/A。今回A,Bの大きな5変数gcdや逆元計算を行わず、保存済みinteger primitive C₁の全係数を使用した。

## balanced degree4 pair

二shape写像は、保存された79項・41項のprimitive分子N₇₉,N₄₁を用いて

\[
 ZN=Z_0\frac{N_{79}}{20F_nJ}
       +Z_0^2\frac{5N_{41}}{4r(r-1)^2J^2}.
\]

s次数はそれぞれN₇₉=3、N₄₁=2、Fn=2、J=1。任意のp(s)についてH_d(p)=A^d p(−B/A)をhomogeneous polynomialとして定義し、

\[
 J_H=H_1(J),\quad F_H=H_2(F_n),\quad
 N_{H1}=H_3(N_{79}),\quad N_{H2}=H_2(N_{41}),
\]

d=r(r−1)²と置く。すると

\[
 \boxed{U=Z_0dJ_HN_{H1}+25Z_0^2F_HN_{H2}},\qquad
 \boxed{V=20dF_HJ_H^2},\qquad ZN=U/V.
\]

U,Vは共に(A,B)に関して4次。実際のsmall79/41係数を使い、元の分子・分母をsについて4次homogenizeしたものがこのU,Vに一致するという**全係数のsymbolic identity**を検証した。大きなC₁係数をsymbolicに代入して展開してはいない。

generic fieldではV≠0。元のdense source chartでFn,J,r,r−1は単元で、Aはbirational linear-C₁回収で非零だからである。このgeneric条件と、後述のspecial algebraでVが単元かどうかは区別する。

## 有限自由モデルと非零determinantの意味

Gのleading coefficientをℓとする。Z_(32003)[R,Z₀,Z₁,Z₂]を下記旧点の係数maximal idealで局所化し、ℓを単元としたringをB_locとする。ℓのresidueが非零なのでこの局所化は適切である。

\[
 \mathcal A=B_{\rm loc}[r]/(G/\ell)
\]

は**明示的なmonic10次多項式**による商で、基底1,r,…,r⁹に関してfree rank10。一般初期idealのbase changeや未証明のflatnessを使っていない。A,B,U,Vはこの有限自由algebraに定義される。C₁の係数は整数で、mod G reductionに使う分母はℓだけ。

10個の列を

\[
 W_j=U^jV^{9-j}\bmod G,\qquad j=0,\ldots,9
\]

とし、この基底での係数を10×10行列Mとする。Mのspecial residue行列のdeterminantが非零なら、generic K上のdeterminantも非零である。generic LではVは非零なので、V⁻⁹を全列に掛けるK-linear automorphismにより、1,ZN,…,ZN⁹も独立。従って

\[
 \boxed{K(ZN)=L}.
\]

**special algebraではAやVを一度も逆にしていない。Vのspecial単元性は仮定も要求もしない。** generic fieldでV≠0というsource条件と、有限自由係数determinantのspecial非零性だけを使う。

## 同じ旧β点と全数値certificate

点は`Q_NEW_FAMILY_OCTIC_BASE_DOMINANCE_20261103_v63`の保存済みindex0であり、元のsource pointは(n₀,n₁,n₂)=(7,9,9)、scale1、base=(−1,1,1,−2)、Sp₄座標は保存値のまま。新しい点を選んでいない。元の16β有理数を全て読み、全denominatorが32003で単元なことを確認してから、同じquadratic C/D係数と4×4 pencil determinantからratioを計算した。

\[
 F_4=30651,\quad(p_0,p_1,p_2,p_3,p_4)
 =(24264,30976,15826,19536,11718)\pmod{32003}.
\]

全quartic trace ratioとF₄²=400p₀+1600p₄の規格化もこの点で確認した。係数baseは

\[
 (R,Z_0,Z_1,Z_2)=(1079,6642,28042,8843).
\]

Gのleading residueは**15049≠0**、monic化後の昇冪係数は

```text
26797,3341,1505,28301,462,12696,30721,19283,29234,5081,1.
```

保存C₁のconstant/linear全10314/10698項をこのbaseにevaluateし、monic Gでrだけをreductionした。U,Vのresidueも全10係数を保存した。行順1,r,…,r⁹、列順UʲV⁹⁻ʲの全100成分を用い、

\[
 \boxed{\det\bar M=21349\pmod{32003}\ne0}
\]

を得た。従ってgeneric QQの係数determinantは厳密に非零で、ZNは10次primitive elementとなる。これは同じold pointにおける**係数modelの非零certificate**であり、source点での零を関係式の恒等性と扱う方法ではない。

## 構成的なrの回収と次のfield-inclusion theorem

generic MについてM·a=e_r、M·b=e_1をCramerの有限なrational circuitで解けば、

\[
 r=\frac{\sum a_jZN^j}{\sum b_jZN^j}
\]

となる。両辺の共通V⁹が消えるからである。sは−B(r)/A(r)、qとBは既存Möbius式で戻る。この存在・構成recipeに巨大なnormの展開は不要だが、generic inverse行列の全symbolic係数は今回展開していない。

全actual NEW target体の指数1への帰結には、別のfield-inclusion論証を用いる。既存の全16座標scale検証はdegree14不変式のsource形をξ⁷×shape、F₄をξ²×shapeとする。全不変式の定義次数下界26は任意の非零Inv₁₄がNEWで消えないことを保証する。shapeを回収した後、その比からξを回収できるという条件は揃うが、**本計算packetの主張はZNのprimitivityまで**とし、そのfull target-field結論は別に明記・監査する。NEWの具体的定義方程式を得たわけではない。

## 再現・保存

```bash
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_ZN_primitive_krylov_20261104_v63.py \
  --output work/ZN-krylov-replay
```

新規outputを指定する。one old point、one prime、60秒/2GB/32MiBの上限で約1.21秒、保存約8.5KB。gateやdetが失敗しても別点・別素数へ自動移動しない実装で、今回は全gateが通過した。

一次保存は`reports/Q_NEW_FAMILY_ZN_PRIMITIVE_KRYLOV_20261104_v63/`。元のβQQ/Fp、C/D、pencil、base ratio、special G、C₁全係数のevaluate値とreduction、U/V中間residue、全100行列成分とdet、入力hashを保存する。source26、new Inv₁₄展開、N norm、target equationの消去、new point/prime probeは未実行。
