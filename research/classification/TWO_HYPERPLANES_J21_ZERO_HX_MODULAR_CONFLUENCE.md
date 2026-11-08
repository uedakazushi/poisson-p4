# Jordan (2,1)、(L,M)=(0,h+x) の u 依存全源を覆う modular 実変形

2026-11-08、PRIMARY `E_reduced_nonnormal`。元の全六源の exact normalizer と、
非零半単純 cofactor を持つ actual polynomial family を使う。新成分は零。
以下は手導出であり、全係数の有限照合と独立 review を別に要する。

## 1. n6 非零の全 chart と分母の無い全六源 normalizer

`TWO_HYPERPLANES_CUBIC_SUPPORT.md` §6.1 および
`CLASSIFICATION_TWO_HYPERPLANES_CUBIC_NILPOTENT_TRANSLATION_COMPONENT_20261107.json`
の全必要十分 source を使う。`S=h+x` と置き、exact cofactor は

```
V0=v∂u+S∂w.
```

元の19係数源の三次式の u 係数は `-n6*S²`。
n6 非零なら `(u,v)->(a*u,a*v)`, `a=1/n6` は V0 を exact に保つ実 PGL である。
vertical determinant が a² なので変更後の F は `a² F(h,x,u/a,v/a,w)`、
従って n6 は a*n6=1 へ移る。唯一反転するのは当該 chart の n6 で、n6=0 全枝は別に保持する。
正規化後の必要十分な六自由変数は `n0,n2,n3,n4,n5,n16`、他の全係数は同節の式(30)で決まる。

以下の全三角 PGL を取る。

```
U=u-n0*h-n3*x,  V=v,
W=w+(n16+n2)*h+(n16+n4)*x+n5*v.                              (1)
```

det=1。逆は `u=U+n0*h+n3*x`, `v=V`,
`w=W-(n16+n2)*h-(n16+n4)*x-n5*V`。source parameter を一つも反転しない。
V0 も exact に保つ。`k=n16`, `K=-S*u+v*w` と書けば全括弧は

```
{x,u}=xv, {x,v}=0, {x,w}=xS,
{u,v}=K-k*S*v,
{u,w}=-K+k*S*w+2k*xv,
{v,w}=xv+k*S².                                               (2)
```

全 k の actual 五次式は `h*x*F0`、

```
F0=-S²u+Svw+xv².
```

これは u について一次で、係数 S² と残部 v(Sw+xv) は互いに素なので既約である。
五つの余分な source 係数を指定点で零に置いた証明ではなく、全逆のある actual normalizer である。

## 2. 二次 base tensor と Hamiltonian/modular field

e を新しい parameter とし `T=S+e*w`, `Ke=-T*u+v*w` を使って、
h-Casimir quadratic tensor Qe を

```
{x,u}=xv,  {x,v}=e*xv,  {x,w}=x*T,
{u,v}=Ke,  {u,w}=-Ke,  {v,w}=xv                              (3)
```

と定義する。これは全 e の実 strict Poisson family である。
同じ係数の全350式も直接照合する。実五次は

```
De=h*x*[xv²+(T+e*v)*Ke].                                    (4)
```

体積形式 `dx∧du∧dv∧dw` について Z_e^i=Σ_j ∂_j Q_e^{ji} と定めると

```
Ve=v∂u+e*v∂v+T∂w=Qe#dlog x,
Ze=-2e*x∂x+(2v-w-e*u)∂u+(e*v-T)∂v+(2h+4x+2e*w)∂w.          (5)
```

Hamiltonian field Ve と modular field Ze は Poisson field である。
`div Ve=2e` が定数なので `[Ve,Ze]=0`。従って

```
Pe,k=Qe+k*Ve∧Ze                                               (6)
```

は全 e,k の実 Poisson tensor である。全係数は QQ[e,k,h,x,u,v,w] の多項式。
e=0 の Pe,k は (2) に正確に一致する。
この modular field は h=1 の実 affine 体積に対するものを斉次化した field であり、
divergence-zero canonical homogeneous lift や全 projective stabilizer を代用していない。

## 3. 一般 e は既受理の非零半単純全支持に入る

`c=1+2ke` と置く。(6) の x 行は厳密に `c*x*Ve` である。
e 非零で fiber の actual affine/PGL 変更

```
R=u-v/e,  A=v,  B=w+(h+x)/e                                 (7)
```

を取ると、Ve は `e*(A∂A+B∂B)` へ送られる。
det=1、逆は `u=R+A/e,v=A,w=B-(h+x)/e`。
従って e*c 非零の開上で非零 tensor scale を除く exact cofactor は diag(0,1,1)、
零固有空間 R への平行移動は零である。

F0 は既約、h*x*F0 は非零なので、各固定 k に対して一般 e でも cubic F_e は既約、rank は四である。
これは proper factorization multiplication map の閉像の補集合の開性を使う。
e=0 で c=1 だから c の零点が原点の近くの全曲線を潰すことはない。
この一般 e の実構造は `TWO_HYPERPLANES_CUBIC_SUPPORT.md` §4.2–4.3 の
全 diag(1,1,0) 支持定理に入り、既知 MPS40 の有限閉和に入る。
Pe,k は e=0 に正則なので、全 k の (2) もその閉和に入る。

ここで rank-one binary form の開だけを抜き出したり、既知一点の帰属を逆向きに使ったりしない。
適用する §4 の定理は全15係数・全binary-rank/complementを含む既受理の全source定理である。
e=0 と c=0 は当該正規化 (7) の補集合だが、端点で必要なのは generic e の帰属と閉性だけであり、
同じ actual family が全 k と全 source 端点を保持する。

## 4. 有限照合と全sourceへの結論

元の必要十分な全六源 binding、normalizer と全逆、全150 tensor 係数、
Qe と Pe,k の strict10/全350 Jacobi、Ve/Ze と Lie/可換恒等式、
cofactor の exact semisimple 変更 (7)、実 Pfaffian と全端点を一つの finite job で照合する。
独立 review が受理されれば、n6 非零の全 source とその全 PGL 飽和閉包が MPS40 に入る。
新しい anchor、局所法線、二次余核の再計算、高次 jet を必要としない。

n6=0 全枝、`(L,M)=(h+x,h)` の p=0、平面 1+1+3 の残三配置、
他の nilpotent source と全 E 完成は別に残す。
