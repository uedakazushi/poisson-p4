# 二conic実衝突の親はそれぞれT16とT17

2026-11-08。全源SNC論証で既知MPS有限閉和へ入れる二つのhelperを、
受理済みcompleted MPS sourceと全括弧で同定する短い補足である。
元のproducer/report/job/manifestは変更しない。

## T16：contact31_chordの親

入力は `reports/MPS_S61_RAW_CATALOGUE_20261103_v62/T16.json` の
全period rho0,rho1,rho2と全edge tau0,tau1を持つhomogeneous bivectorである。
原座標を(x0,x1,x2,x3,x4)=(A,B,C,D,H)と書く。
EをEuler場とし、

```
L_A=L_H=0,
L_B=(rho1−rho2)B/2,
L_C=(rho1+rho2)C/4,
L_D=−(rho1+rho2)D/4
```

とする。これは元のH行をHで正確に割ったlinear fieldである。
P_raw−E∧LはH-Casimirとなり、その全六括弧は

```
{A,B}=p AB+tau0 H², {C,D}=q CD+tau1 H²,
{A,C}=k AC, {A,D}=−k AD,
{B,C}=−k BC, {B,D}=k BD,
p=(rho2−rho1)/2,
q=rho0+(rho1+rho2)/2, k=−(rho1+rho2)/4.
```

従ってchord helperへのexact parameter bindingは

```
(rho0,rho1,rho2,tau0,tau1)
=(q+2k,−p−2k,p−2k,p u,q v).
```

Euler gaugeは同じprojective bivectorであり、座標順の変更も固定PGLである。
全五parameterの一般点でこの式は受理済みT16 sourceへ入るので、
閉性からhelper全体がT16に入る。
`TWO_QUADRICS_CHORD_SNC_CONFLUENCE.md` の全parameter実族と全六端点により、
**contact31_chordの全六磁気源と全PGL飽和はT16に入る**。

## T17：contact31_tangent_at_simple_intersectionの親

入力は `reports/MPS_SMOOTHING_COMPLETION_20261103_v62/T17.json` である。
未完成の単純edge和ではなく、tau0 tau1 H²/rho1項を含むcompleted sourceを使う。
同じ原座標順において

```
L_A=L_H=0, L_B=−rho1 B,
L_C=−(rho1−rho2)C/2,
L_D=−(rho1+rho2)D/2
```

を引くEuler gaugeにより、全六括弧は厳密に

```
{A,B}=F, F=p AB+uH²,
{A,C}=a AC, {A,D}=b AD,
{B,C}=−a BC, {B,D}=−b BD,
{C,D}=q CD+vF,
p=rho1, a=(rho1−rho2)/2, b=(rho1+rho2)/2,
q=rho0+rho2, u=tau0, v=tau1/rho1.
```

ここではa+b=pとvFの積関係がそのまま成立する。
coupled helperからの逆bindingは

```
(rho0,rho1,rho2,tau0,tau1)
=(q−p+2a,p,p−2a,u,pv).
```

p≠0の稠密開で全completed sourceと一致する。
helper自体はp=0にも多項式的に定義されるので、その補集合はT17の閉性で含まれる。
実際のtan-simple衝突ではp=1を固定しており、原source分母rho1は全epsilonで単元である。
従って
`TWO_QUADRICS_TAN_SIMPLE_SNC_CONFLUENCE.md` の全実族と全六端点から、
**contact31_tangent_at_simple_intersectionの全六磁気源と全PGL飽和はT17に入る**。

## 依存と範囲

各completed sourceの閉包がnamed MPS componentであることは、
`research/notes/MPS_S61_CATALOGUE_20261103_v62.md` §§1–2の受理済みsource認識への依存である。
今回の新しい部分は上の全period/edgeの六括弧と明示Euler fieldの直接比較だけである。
一点のPfaffian一致から成分所属を推測していない。
全zero磁気境界を含む実衝突は既存の別producerと独立reviewで検証される。
このnamed同定が無くても、両SNC helperの既約性と全SNC成分定理に基づく元の被覆証明は独立に有効である。
