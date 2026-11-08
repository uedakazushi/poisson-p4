# Smooth cubic scroll rank零：全十二 Jacobi source の四閉源への純還元

2026-11-08、PRIMARY `E_reduced_nonnormal`。被約 support の全 case split。
入力は accepted abelian-conormal 十二源と元350行の独立二十四 generator。
新 GB、native、rank、minor 又は point sample は使わない。
この段階では実 D の因子、rank四存在、既 component 帰属を主張しない。

## 1. 記号と二十四元式

`CLASSIFICATION_CUBIC_SCROLL_ABELIAN_DIMENSION_20261107.json` cases[0].generators をそのまま使う。
番号 g1,…,g24 は保存順である。次と置く。

```
p=r0, a=r1, s=r2, q=r3, v=r4, f=r5,
r=r6, h=r7, u=r8, w=r9, b=r10, t=r11,
δ=b−a, η=t−2s, ell=q+h.
```

旧 full350 の quadratic row rank は二十四で、保存generatorはその元行から選んだ basis。
従って二十四式の零と元全350式の零は全 source 上で同値である。
以下は reduced point support の議論であり、scheme ideal の radical性を仮定しない。

## 2. 六次元 linear source L

まず元 g19=(t−s)(t−2s)、g1=ηf、g9=ηu、g20=ηδ を使う。
η非零では t=s、f=u=δ=0、s非零。
g24 から h=0、g3 から w=−v。従って

`L: t=s, b=a, w=−v, f=u=h=0`

へ入る。自由変数は p,a,s,q,v,r の全六つ。

η=0 かつ v非零では、g4=v(−3δ+4f)、g8=v(δ−f) から f=δ=0。
g5 より w=−v、g3 より s=0、g13 より h=0、g12 より6vu=0。
この whole complement も L に入る。
L 自身へ全二十四式を代入すると全て零である。

## 3. 残余は一つの九式へ正確に戻る

残る η=v=0 では t=2s、v=0。
g2−g3=(δ+2f)(a−2f)、g6+6g3=−2(a−2f)(δ+f)。
よって (a−2f)f=(a−2f)δ=0。
さらに g22 から s(a−2f)=0、g7 から w(a−2f)=0。
a−2f非零では f=δ=s=w=0、g24 より u=0、g12 より h=0で、再び L。

したがって L 外では

`t=2s, v=0, a=2f, b=2f+δ`

としてよい。元二十四式の零は、この上で正確に次の九式の零と同値である。

```
R1=s*w−f*δ,
R2=f*u,
R3=s*h−δ*u,
R4=p*f+s*(q+h),
R5=u*w−f*h,
R6=p*w+δ*q+(δ−2f)*h,
R7=p*u+s*r−2u²,
R8=−p*h−(δ+2f)*r+2q*u+4h*u,
R9=h*(q+h)−r*w.
```

例えば g10=2R4−4R2+R3、g11=−2R4+2R2−R3、
g12=−2R6+2R5、g14=R6−R5、g15=R6+R5、
g16=R7、g17=R8、g18=R9、g24=R3。
g2=g3=R1、g6=−6R1、g21=6R1、g23=−4R1、他は零。
従って逆方向にも定数二／三だけで R1,…,R9 を元generatorから回収出来る。
source parameter を反転する ideal 同値ではない。

## 4. 九式の全 unit／zero 補集合

f非零では R2,R5から u=h=0。
R7,R8,R9から sr=(δ+2f)r=wr=0。
もし r非零なら s=w=0、δ=−2fとなり、R1が2f²=0を与えて矛盾。
従って r=0で、残式はちょうど

`rank [[f,s],[-q,p],[w,δ]] ≤1`。

この全閉 Segre coneを F と書く。rank一 matrixの全外積表示があるので irreducible、dimension四。
f非零開は dense、zero complementも同じ閉源に保持する。

f=0,u非零では R5よりw=0。R3,R7をR8へ戻すと2u(q+h)=0で ell=0。
残式はちょうど

`rank [[u,s],[-r,p−2u],[h,δ]] ≤1`。

この全閉 Segre coneを U と書く。ここもdimension四で u非零開の閉包全体を保持する。

f=u=0,s非零では R1,R3,R4,R7から w=h=ell=r=0。
残 p,δ,s は自由。この全枝は F と U の両方の閉源へ入る。

f=u=s=0では残式は

`p*w+δ*ell=0, p*h+δ*r=0, h*ell−r*w=0`。

すなわち

`rank [[h,w,−δ],[r,ell,p]] ≤1`。

この全閉 Segre coneを S と書く。dimension四で全 rank零 matrixも含む。
以上は反転した f,u,s,η,v,a−2f の全零補集合を尽くす。
逆に F/U/S の各全 rank≤1 source で R1,…,R9 は恒等零。
Lと合わせ **全被約 Poisson source は L∪F∪U∪S** である。

## 5. 全 source を覆う具体的 outer-product map

順序は元の(r0,…,r11)。Lは全六 parameter (P,A,S,Q,V,R) で

`(P,A,S,Q,V,0,R,0,0,−V,A,S)`。

各 Segre cone は全五 parameter (A,B,C,D,E) の outer productから次で覆う。

```
F: (B*E, 2*A*D, A*E, −B*D, 0, A*D,
    0, 0, 0, C*D, 2*A*D+C*E, 2*A*E),

U: (B*E+2*A*D, 0, A*E, −C*D, 0, 0,
    −B*D, C*D, A*D, 0, C*E, 2*A*E),

S: (B*E, 0, 0, B*D−A*C, 0, 0,
    B*C, A*C, 0, A*D, −A*E, 0).
```

F/U は3×2 matrix=(A,B,C)^t(D,E)、Sは2×3 matrix=(A,B)^t(C,D,E)。
C上の任意の rank≤1 matrix はこのように因数分解出来、零matrixも含む。
chosen chart や generic nonzero coordinatesだけの mapではなく、全 source の全 parameter像である。
余分なC* redundancyは source を失わせず、有限全係数評価に問題を生じない。

次の有限義務は、既 complete30×12 source matrix からこれら全 map の actual tensorを作り、
元 full350 Jacobi と全126 actualDを一度だけ全 parameterで評価することである。
最弱目標は全 actual factor／rank四空性／essential-three 又は既知 whole family。
元原式の compatibility を省いた点の集まり、factorized symbolからの推測、GB／minor追加へは進まない。
