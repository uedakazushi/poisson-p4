# Stable nil22 二全源は T21 の全 factor-flag 閉包の外

2026-11-08。PRIMARY `E_multiple_hyperplane`。新計算を行わない手証明。
全六磁気係数、零磁気、反復する二つの marked hyperplanes を含む。
本稿は独立読解待ちであり、他の親の不存在又は E 全完了を主張しない。

## 1. 元の全源と共通の二つの場

座標を `(h,q1,q2,x,y)`、`p=(x,y)` とする。
`NIL22_STABLE_PROJECTIVE_ACTION.md` と `NIL22_IMAGE_SOURCE_REVIEW.md` の
全源は h-Casimir homogeneous lift で

```
P#dh=0, {qi,pj}=Mij, {x,y}=0,
{q1,q2}=R+ h*(V2*q1-V1*q2),
R=m0*h²+m1*hx+m2*hy+m3*x²+m4*xy+m5*y².
```

`M` は base 変数 `(h,x,y)` だけを含み、両 stable 源とも

```
M|h=0 = -p*p^t/3.
```

従って column divergence `Z_h,j=sum_i partial_i Pij` は

```
Z_h|h=0 = V := x*partial_q1+y*partial_q2.                 (1)
```

実際、q 成分では q1q2 項の微分は h の倍、M の一行の二次部分の発散は
`-p_i`、他の成分は零である。R は q を含まないので (1) に寄与しない。

二源の実 D はそれぞれ

```
line_conic: h²*L*Q,  L=x-y-12h,
 Q=36h²+3h(x-y)+2xy   (rank3),
three_lines: h²*(x+3h)*(y-3h)*(x-y-3h).
```

全一次因子は base の一次式である。L が h と比例しない任意の一次因子なら、
全体に正則な linear log-Hamiltonian を

```
X_L=P#dlog(L/h)=P#dL/L
```

とする。因子の Poisson 接性と h-Casimir 性から分子は L で正確に割れる。
`L!=h` の h 上の稠密開で M の上式を割れば

```
X_L|h=0 = -V/3.                                         (2)
```

両辺は正則な linear field の制限なので、L=0 と交わる点にも (2) が延長する。
`X_h=0` と定める。磁気係数、L の係数、又は座標 x,y を反転していない。

## 2. 二超平面が合体しても定義される proper flag 条件

Poisson tensor と Poisson 超平面 H に対する H-Casimir homogeneous lift は一意である。
差が `E wedge Y` で H 行が零なら `H*Y-Y(H)*E=0` より差は零だからである。
H のある係数が非零の chart で H を第一座標に取れば、Euler gauge とこの lift は
有限の線形係数抽出・正確な一次式除算で正則に得られる。
発散を `Z_H` と書く。

もう一つの非零 Poisson 一次式 W を取る。**H,W が同じでもよい。**
`P_H#dW` は W の倍である。H,W 独立の場合は W 上で projective 接性と
H-Casimir 性を合わせて証明でき、W=H の場合は分子自体が零である。
一次式 W を掛ける `degree1 -> degree2` の写像は全 W!=0 で一定階数なので、
その像は vector subbundle となる。従ってその一意な商

```
X_{W/H}=P_H#dW/W
```

も W の係数が非零の chart で正則である。これは W!=0 の affine 空間へ
局所化するという意味ではなく、全係数の正確な割り算である。
特に W=H では `X_{H/H}=0` で、二 marked factors の衝突に pole はない。
この議論は reduced incidence/support 上で十分であり、source の非被約構造を仮定しない。

以下を ordered factor flag の閉条件とする。

```
Z_H+X_{W/H} is projectively zero on W.                    (C)
```

これは `E wedge (Z_H+X)` の全係数が W で割れるという多項式条件である。
上の正則性から、独立 flag のみでなく H=W 上でも閉条件である。

tensor と因子の graph を
`P(H) x P(W) x P(cubic)` の proper 積の中で閉じる。
乗法 `H*W*cubic=D` は proper であり、D!=0 の任意の boundary tensor は
極限 factor flag を持つ。各 flag 上で (C) は延長する。
従って D の cubic 因子が可約になる場合、又は H=W となる場合も、
その全ての因子割当てを調べればよい。

## 3. T21 の完成全源は条件 (C) を満たす

既独立受理済み `TWO_LINES_CUSP_TANGENT_SECANT_GLOBAL_PARENT_SEPARATE_REVIEW.md`
§4 の全五 parameter 係数恒等式を再利用する。
完成 T21 は H-Casimir gauge で

```
P_H=Jac(f_ABC)+p0*W*partial_W wedge E_ABC,
f=(p1*ABC+t0*C³-t1*B³+t2*A³)/3,
Z_H=p0*E_ABC-3*p0*W*partial_W,
X_{W/H}=-p0*E_ABC.
```

従って `Z_H+X=-3*p0*W*partial_W` は W 上で零である。
period/edge を一つも割らず、全 source 上の literal identity である。
既受理の完成 source 定理から PGL 飽和閉包は T21 全成分を与える。
従って同全成分の全 rank4 boundary は §2 の何らかの ordered factor flag で
条件 (C) を満たす。

## 4. stable target の全 factor flags が (C) に反する

旧 h-Casimir 代表から H-Casimir 代表への変換は

```
P_H=P+E wedge X_H,
Z_H=Z_h+5*X_H modulo Euler.
```

後式は `delta(E wedge X)=5X-(div X)E` の直接微分である。
また `X_{W/H}=X_W-X_H` なので、条件 (C) の場は

```
Z_H+X_{W/H} = Z_h+4*X_H+X_W modulo Euler.                (3)
```

target D の全一次因子は §1 の base 一次式である。
極限 flag の H,W はそれぞれ h 又は h と比例しないその一次因子である。
(1),(2) を (3) に入れると、旧 h 上の場は次の非零倍の V になる。

| H | W | (3) on h=0 |
| --- | --- | --- |
| h | h | V |
| h | non-h factor | 2V/3 |
| non-h factor | h | -V/3 |
| non-h factor | non-h factor | -2V/3 |

実 D の重複度から不可能な組も表に残した。従って factor の選択漏れはない。
特に `(H,W)=(h,h)` だけへ衝突を限定していない。

W=h の場合は h 上で、W!=h の場合は `h=W=0` 上で `(x,y)!=0` の点を取る。
後者の base 直線条件は x,y の非零解を持ち、q1,q2 は自由である。
V は base 成分が零で、Euler は非零の base 成分を持つから、この点で
V は Euler と比例しない。従って表の各非零倍も W 上 projectively 零ではない。
条件 (C) は全 flag で失敗する。

## 5. 結論と残余範囲

`stable_line_conic` と `stable_three_lines` の**全六磁気源**は T21 の全成分の外にある。
磁気の正規化・非零不変量の仮定・source 一般点の選択は不要であり、
零磁気と全磁気分母補集合を含む。

この whole T21 除外は、旧 uniform collision の特定 frame 不整合からは導いていない。
衝突 factor flag 上の正則な intrinsic field 条件を用いた別の全源論証である。
T36/T37 への陽性所属、全既知親の不存在、新成分性、又は E 完了は結論しない。
既保存 T38 Chow 条件と合わせても、その scope 以外へ外挿しない。

新 CAS、rank、minor、operator、point、normal jet は実行していない。
本稿と新 manifest のみを作成し、既 producer/report/job/独立 review は全て不変である。
