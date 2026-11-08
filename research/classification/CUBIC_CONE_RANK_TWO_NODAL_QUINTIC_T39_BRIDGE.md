# 残cubic-cone一軌道とnodal elliptic quintic：T39への最弱線形認識

2026-11-08、PRIMARY E。新local obstruction、jet又はpoint列は使わない。
本稿は次の**一つの有限条件付き認識**を証明する。

canonical代表の既知quinticをFとし、全150二次係数に対する

`L_F={P : div(P)=0, P*grad(F)=0}`

が一次元なら、この代表とそのwhole sourceは閉成分T39に入る。
必要なのはこの一つの全線形kernelと以下のliteral恒等式だけである。
条件を満たすかは同名finite reportで判定し、次元を仮定しない。
全pointの新componenthoodを主張するものではない。

## 1. 前の有限義務を正確に停止する

同じ代表についてD100、Euler24、actual M49が得られた。
真正150直和M49+R+U100のdetは6103515625000000000000。
一本normalのquadratic Jacobi余核は全350成分が零で、[D|JacR]はrank100だった。
従ってこの有限testはcomponenthoodを閉じず、契約通りcubic／高次jetへ進んでいない。
`CLASSIFICATION_CUBIC_CONE_RANK_TWO_ONE_NORMAL_QUADRATIC_20261108.json`
を零という結果のまま保存する。

## 2. 全zero曲線のexplicit normalization

前noteでSing(F)_red=S、c-vertex外のS上restrictionは
`5(1−t^5+t*eta) ∂t wedge ∂eta`と得た。
その零曲線のclosureは

`E0=[s^4*t:s^3*t²:s²*t³:s*t^4:t^5−s^5]`。

五sectionは独立で共通零を持たない。st≠0ではx1/x0=t/sが逆なのでbirational。
二端点だけが同じc-vertexへ写る。c=1の局所座標で一方の接線はx0軸、他方はx3軸。
従ってordinary nodeが一つだけで、他の点はratioによる局所逆からsmooth。
像はintegral degree五、算術種数一、全P4をspanするnodal elliptic-normal quinticである。

この曲線上でvanishする五つの独立quadricは

```
q0=x0*x2−x1²,
q1=x0*x3−x1*x2,
q2=x1*x3−x2²,
r1=c*x1−x2*x3+x0²,
r2=c*x2−x3²+x0*x1.
```

その五行Jacobianのdeterminantは正確に

`det Jac(q0,q1,q2,r1,r2)=−2F`。

Fは前のwhole-source noteの
`c*q0*q2+q2*(x2*q2−x3*q1)+q0*(x0*q1−x1*q0)`。
finite認識ではこの全polynomial恒等式も直接確認する。

## 3. regularなWeierstrass smoothingを明示する

Weierstrass nodal cubicを `y²=x³+x²` とする。
u=(t+s)/(t−s)、x=u²−1、y=u(u²−1)が正規化で、
smoothなpoint at infinityはt=s。L=O(5*infinity)の基底は(1,x,y,x²,xy)。
上のP1上のsectionに戻すと、次の**constant可逆linear map**となる。

```
z0=c+5*x0−10*x1+10*x2−5*x3,
z1=−4*x0+12*x1−12*x2+4*x3,
z2=4*x0−4*x1−4*x2+4*x3,
z3=−16*x1+16*x2,
z4=16*x1+16*x2.
```

実際z0=(t−s)^5、残四つはz0*(x,y,x²,xy)である。
正規化chartの分母をactual PGL mapと取り違えず、最後に使うmapは上のconstant GL5のみ。

parameter eについてsmooth cubic族

`y²=x³+x²+e`

を取り、同じinfinityとLを使う。discriminantは非零定数倍の
`e*(4+27e)`なのでgeneric eでsmooth。degree五完全linear systemの像Eeを
P4_zで表す五つの独立quadricsは

```
g0=z0*z3−z1²,
g1=z0*z4−z1*z2,
g2=z1*z4−z2*z3,
g3=z2²−z1*z3−z0*z3−e*z0²,
g4=z2*z4−z3²−z1*z3−e*z0*z1.
```

z0≠0でz3=x²,z4=xyと元cubicへ戻る。
z0=0には唯一のinfinity点しかなく、全degree五complete-series modelである。
特殊e=0の五quadric spaceは上のconstant mapの下で(q0,q1,q2,r1,r2)のspaceと一致する。
これは全五列のconstant invertible basis changeとして有限receiptへ保存する。

## 4. smooth側secant quinticの一次資料とspecializationの範囲

一次資料はT. A. Fisher, *Invariants for the elliptic normal quintic*,
https://www.dpmms.cam.ac.uk/~taf1000/papers/invenq.pdf 、§9、印刷33–34頁。
2026-11-08に全文PDFを取得し、§9–10を直接読解した。
PDFのSHA256は有限入力receiptへ固定する。
同§9の式(34)とその前のsecant constructionは、smooth elliptic normal quinticの
五quadric basisのJacobian determinantがsecant hypersurfaceの式であることを与える。
quadric basis changeはdetに非零constantを掛けるだけ。

ここではsmooth generic Eeにだけこの命題を使う。
nodal special curveへsmooth-only theoremを直接流用していない。
F_e=(det Jac(g0,...,g4))を上のconstant GL5でxへ戻すと、F_eはeのpolynomial族で、
F_0は前節のnonzero scalar倍のFである。
従ってspecialization後も実quinticは零にならない。

generic Eeのsecant quinticは既存normal/q5,1 theoremのT39 actual退化因子である。
同smooth側のq5,1構造をπeとする。canonical strict lift Peは
`div(Pe)=0, Pe*grad(F_e)=0`を満たす。
存在とT39帰属は既受理のnormal actual divisor theorem／elliptic secant構成だけを使い、
未知のspecial zero-scheme分類を仮定しない。

## 5. proper incidenceと一次元kernelによる全tensor認識

全canonical quadratic liftのprojective spaceとe-line上で

`I={(e,[P]): div(P)=0, P*grad(F_e)=0}`

はclosed incidenceで、e-lineへの射影はproperである。
generic側の([Pe],e)のclosureはspecial fiberに非零projective point[P0]を持つ。
PeはT39内なのでP0も閉T39内にあり、L_Fに属する。

もし**元全150係数のL_Fが一次元**なら、既知の我々のcanonical代表がその非零generatorなので、
P0は正確にそのprojective classである。
これがT39への陽性包含を証明する。Fの形だけからtensorの帰属を推測したのではない。
specialized covariantが非零だと仮定する必要もなく、projective propernessが非零limitを保証する。

前のwhole-source proofでλ3/5のsemisimple unit全源はこの代表の一actual軌道、
zero二cohom枝とnil全六fiberはその軌道閉包に入る。
したがってこの陽性認識が通れば**λ3/5の全sourceはT39**へ戻る。
λ2/5全sourceは既T38、残familiesは既約Dで空という前の全補集合も保持する。

L_F次元が一を超えればこの認識は未閉と保存し、同じoperatorの追加anchor列又は
局所高次計算へ自動継続しない。全global component数の増分は零のままである。
