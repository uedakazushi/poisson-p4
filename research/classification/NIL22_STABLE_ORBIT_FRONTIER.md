# stable 二源：全 toric 親と新103親を除く構造的境界

2026-11-07、`E_multiple_hyperplane`。このノートは既存全源・全行列証明書を
読み直した上での数学的帰結であり、新しい source 所属を証明したものではない。
固定 confluence の不整合を一般化せず、群作用の閉条件と proper な flag incidence を使う。

## 1. 正確な軌道次元と、全 T01 に対する障害

`CLASSIFICATION_NIL22_IMAGE_STABLE_CUBIC_20261107.json` の二つの固定 witness は、
全150係数上の Euler 鉛直・GL5・scale 行列について次の正確な階数を持つ。

| 全六磁気源 | 鉛直階数 | 全作用階数 | projective PGL 軌道次元 | source 像下界 |
|---|---:|---:|---:|---:|
| stable_line_conic | 24 | 47 | 22 | 23 |
| stable_three_lines | 24 | 46 | 21 | 23 |

projective 軌道次元は全作用階数から `24+1` を引いた値である。
three_lines を orbit22 と読み違えてはならない。両源は polynomial な全六磁気源であり、
各 PGL 飽和像閉包は既約である。表示階数は report の非零小行列による既済の入力である。
このノートで小行列を新たに再計算したとはしない。

MPS の T01 の稠密全源は、独立な五つの一次式を座標とした全 logcanonical tensor である。
任意の period で四次元対角 torus が projective stabilizer に含まれるので、
全 source で projective 軌道次元は高々 `24-4=20` である。
projective Poisson 空間上の PGL infinitesimal action は代数的な vector-bundle map
であり、その階数が20以下という条件は全小行列で定まる閉条件である。
従って **T01 の全閉包**で同じ上界が成り立つ。

上記二 witness の軌道次元はどちらも20より大きい。よって両 witness は T01 に入らず、
二つの全 source のどちらも T01 に含まれない。表示した非零小行列の非消失開集合上でも
同じ非包含が成り立つ。この結果は磁気特殊境界の全点非包含までは主張しない。
別の toric moving frame、period 次数、特殊 profile を選んでもこの全親族障害は変わらない。

## 2. 新103の marked tacnode 条件も全二源を覆わない

既知成分 `R4_GA9_CUSPIDAL_PLANE_QUARTIC_24` の一般の実五次は
`D=H*F` であり、F は A3 tacnode と A2 cusp を持つ平面 quartic、H は
tacnode の二枝に共通の接線である。
この正確な入力は `GA9_CUSPIDAL_QUARTIC_SEPARATE_REVIEW.md` による。

一般 marked pair には `p∈H` があり、p を中心とする正則な平面局所座標 `(H,ell)` で

```
F ∈ (H², H*ell², ell⁴).
```

特に p で一次 jet は零、二次接錐は H² の倍数である。
これらの jet 条件は平面上の点・直線 flag の projective incidence における閉線形条件である。
flag の parameter 空間が proper なので、この存在条件は marked pair の閉包に残る。
以下ではさらに弱い「p は F の特異点で、二次接錐は H² の倍数」を主に用いる。

退化因子の本質変数数が generic と special の双方で3であれば、方向微分行列の kernel は
DVR 上の自由直和となり、正則な基底で共通の平面商へ移れる。
これは `GA9_LINEAR_W_MARKED_PAIR_OBSTRUCTION.md` の冒頭と同じ一定階数の議論である。
従って P4 の PGL 飽和閉包についても、次の平面 flag 議論を適用できる。
generic の一次因子 H は積写像の proper 性と有限 base change により special に延びる。
special D の単純因子へ移る可能性も次に全て含める。

stable_line_conic の実因子を

```
D=H0²*L*Q,
H0=h, L=x-y-12h,
Q=36h²+3h(x-y)+2xy
```

とする。Q は滑らかな conic、H0 は割線、L は Q の接線であり、
`H0∩L∉Q`。H0 と C3=LQ の交点は相異なる三点である。

tracked generic H の special が H0 なら `F0=H0*C3` である。
flag 点 p で C3(p) が非零なら F0 は滑らかで不可能。
零なら C3|H0 は単根なので、F0 の二次接錐には非零の H0*ell 項があり、
H0² の倍数にはならない。従ってこの場合は不可能である。

tracked H が L なら `F0=H0²*Q`。L 上で F0 が特異となり得る点は H0∩L だけで、
そこでは Q!=0、二次接錐は H0² であって L² の倍数ではない。
L∩Q の接点では H0!=0 なので F0 自体は滑らかである。これで全一次因子を尽くした。

stable_three_lines は `D=H0²*L1*L2*L3` であり、

```
H0=h, L1=x+3h, L2=y-3h, L3=x-y-3h
```

の四直線は三本ずつ独立である。
tracked H=H0 なら residual cubic の H0 上三単根により前と同じ矛盾。
tracked H=Li なら `F0=H0²*Lj*Lk` であり、Li 上で F0 が特異となるのは
H0∩Li だけで、その二次接錐 H0² は Li² と異なる。
他の交点では F0 は滑らかである。三本が交わらない条件で全点が尽くされる。

従って、実退化因子の幾何だけで、この二つの全磁気源は新103成分の外にある。
double_A では residual cubic の H0 上二重根があるため、この障害は成立せず、
実際に新103の正則 a=1/4 端点へ入る。別ノートの正の結果と矛盾しない。

## 3. 証明状態と未解決

新しい反復計算はしていない。全作用階数は元の全小行列証明書に依存し、
T01 の稠密 source 認識と新103の一般 marked pair も受理済み入力である。
上の閉条件・proper flag・全 factor limit の論証は今回手で確認したもの。
別作者の数学レビューは未実施である。

両 source の E 帰属は未解決である。T39/T40 の全磁気 graph 診断、
CPLUS_EXT/従属三H+Q の閉境界障害はそれぞれ既存ノートへ依存する。
他の MPS24、未知成分、他の実変形を排除していない。
本結果を新成分性や全 E 被覆として数えない。
