# double_A の全六磁気源は既知平面 quartic24 成分に入る

2026-11-07、PRIMARY `E_multiple_hyperplane`。
全係数の実 Poisson 族と端点恒等式を今回直接確認した。別作者 review は別ノートで記録する。
新しい成分は加えない。既知成分は `R4_GA9_CUSPIDAL_PLANE_QUARTIC_24`。

## 1. 既知の全源と正則な構造境界

入力は `GA9_CUSPIDAL_QUARTIC_COMPONENT.md` と独立
`GA9_CUSPIDAL_QUARTIC_SEPARATE_REVIEW.md` による、全六磁気係数を持つ実源

```
F_a(h,y,w)=(hw+6y²)²−144y³(h+ay),
P_a(R)=partial_x∧U_a+partial_z∧V_a
       +(alpha_a*x+beta_a*z+R(h,y,w))partial_x∧partial_z,
R=r0*h²+r1*hy+r2*hw+r3*y²+r4*yw+r5*w².
```

全七 source の PGL 飽和閉包が上記24次元全空間成分となることは既済の入力である。
公開 `CLASSIFICATION_GA9_CUSPIDAL_QUARTIC_COMPONENT_20261107.json` の
`actual_family_power3` に、U,V,alpha,beta の全有理式が保存されている。
この入力の bytes や旧証明書は変更しない。

`a=1/4` では四次因子が

```
F_(1/4)=h*C3,   C3=h*w²+12w*y²−144y³
```

となる。quartic の既約性の開集合の境界だが、Saito frame の分母 `18a−7` は
`−5/2!=0` であり、実二ベクトルの全係数は正則である。
従って任意の r に対して、この端点自体が同じ既知閉成分に入る。
退化因子の一致だけで帰属を推論しないよう、以下で全150係数を合わせる。

## 2. 固定した実 PGL 変換

target は `CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json` の
`cases[0].full_six_magnetic_source`、座標 `(h,q1,q2,p1,p2)` である。
元 report の `(x1,x2,x3,x4)` を `(q1,q2,p1,p2)` と読み替える。
三次因子は

```
C=(384p1²p2+512hp1²+384hp1p2−9hp2²−432h³)/384,
actual_D=-h²*C.
```

parent base を target base の式で

```
y=−4p1−3h,
w=−48p1−(9/2)p2−54h
```

と置くと、直接 `C3=-864*C` である。二変数線形部の行列式は18。
fiber は parent `(x,z)` から target `(q1,q2)` へ

```
(q1,q2)^T=A*(x,z)^T,
A=[[-1/448,-5/1344],[11/252,-5/252]],   det A=5/24192.
```

全五座標の parent→target 行列は

```
G=[[1,0,0,0,0],
   [0,-1/448,0,-5/1344,0],
   [0,11/252,0,-5/252,0],
   [-3/4,0,-1/4,0,0],
   [-4,0,8/3,0,-2/9]],
det G=-5/435456.
```

従って真の可逆 PGL 変換であり、パラメータに依存する極を持たない。
四本の anchor 括弧の全係数を合わせる定数線形系から A を求め、その全係数一致を確認した。
commutator の fiber 一次項も、この G の全 pushforward に含めて確認する。

## 3. 全六磁気係数の可逆な一致

target magnetic polynomial を
`M=m0*h²+m1*hp1+m2*hp2+m3*p1²+m4*p1p2+m5*p2²` とする。
端点での磁気変換は正確に

```
M(h,p1,p2)=(5/24192)*R(h,-4p1-3h,-48p1-(9/2)p2-54h).
```

これは全二次形式の六次元空間上の可逆線形変換である。
その六行列と逆変換を report に全て保存した。行列式は

```
(5/24192)^6 * 18^4 = 15625/1909592158002053382144 != 0.
```

どの m も反転せず、零係数や零磁気項の補集合も同じ恒等式で含まれる。
この逆変換で r を m の式に置き、`a=1/4+epsilon` とする。
得られる真の実族を

```
Pi_epsilon(m)=G_* P_(1/4+epsilon)(R(m))
```

と書けば、公開 report の target 全源に対して

```
Pi_0(m)=P_double_A(m0,...,m5)
```

が全150係数で恒等的に成り立つ。全六 magnetic と epsilon の積も元の source 式のままであり、
Taylor 係数を独立未知数へ置き換えたものではない。

## 4. 全 Jacobi・実 Pfaffian と閉性

上の Pi の全350 projective Jacobi、さらに全10 strict homogeneous Jacobi は
`QQ(epsilon)[m0,...,m5,h,q1,q2,p1,p2]` 上で恒等的に零。
全係数の分母は `192*(36epsilon−5)` の因子だけで、epsilon=0 では単元である。
affine Pfaffian を直接計算して得られる実五次は

```
D_epsilon=5h/(36epsilon−5) * [h*C+(epsilon/6)*(4p1+3h)^4].
```

特に `D_0=-h²*C!=0` であり、rank4 の target を実際に保持する。
Pfaffian の形式因子だけを比較した計算ではない。

generic epsilon は既済の全七 source に属する。よって任意の m に対して
その実族は既知閉成分内にあり、正則な端点も同成分内にある。
既に端点で全六 source を直接含むので、単一 canonical orbit の稠密性だけを用いる必要もない。
PGL 不変性と閉性により、全 double_A 源の PGL 飽和閉包もこの既知24成分に含まれる。
既知境界点の帰属から別の source の帰属を逆推論した議論ではない。

## 5. 検証と残る範囲

```
python3 scripts/classification.py gate research/classification/jobs/nil22_double_a_cuspidal_limit.json
timeout 90s prlimit --as=1610612736 -- .venv/bin/python code/classification_nil22_double_a_cuspidal_limit.py
```

現行 report は `CLASSIFICATION_NIL22_DOUBLE_A_CUSPIDAL_LIMIT_20261107.json`。
全恒等式は2.002秒、最大 RSS 56,248 KiBで直接確認した。
元の既知成分性と target の全十一源入口は既済の独立証明への依存であり、今回再証明していない。
この新しい全源一致は別作者確認前の producer 証明である。gate 自体は数学を認証しない。

独立受理後、nil22 image の全十一源に残る E 未帰属は stable 二源だけとなる。
それらや nil31-short r=2/3、他の未処理 E/G 支持を解決したとは主張しない。
固定 SS12 源の全非包含診断とは別の数学結果であり、その診断を本証明の依存にしない。
