# Smooth cubic scroll rank零：修正された全被約支持 L6＋Segre五

2026-11-08、PRIMARY `E_reduced_nonnormal`。純証明 V2。
旧四源 cover の誤りは `CUBIC_SCROLL_ZERO_FOUR_COVER_ERRATUM.md` に明示し、旧bytesは保存する。
本稿は正しい全二十四元式から再導出し、旧 fu=0 を使わない。
この段階では新 Segre五の actualD因子／既約quintic空性は未検証である。

## 1. 変数と旧還元の有効部分

```
p=r0, a=r1, s=r2, q=r3, v=r4, f=r5,
r=r6, h=r7, u=r8, w=r9, b=r10, t=r11,
δ=b−a, η=t−2s.
```

元二十四generatorは frozen DIMENSION report の保存順 g1,…,g24。
η非零では g19、g1、g9、g20、g24、g3から順に
`t=s,f=u=δ=h=0,w=−v`。sはこの枝で非零。
η=0,v非零では g4/g8から f=δ=0、g5からw=−v、g3からs=0、g13からh=0、g12からu=0。
両枝は全六linear source

`L6: t=s,b=a,w=−v,f=u=h=0`

に入る。
残η=v=0では
`g2−g3=(δ+2f)(a−2f)`、`g6+6g3=−2(a−2f)(δ+f)`。
a−2f非零なら f=δ=0、g22/g7からs=w=0、g24からu=0、g12からh=0で、再びL6。
これらの推論は旧誤g11を一つも用いない。

従って L6 外では t=2s,v=0,a=2f,b=2f+δ。
以下ではこの全subspaceを保持し、元の fu非零を捨てない。

## 2. 正しい全八式と元24との定数両方向

`P=p−2u, L=q+h` と可逆 linear change を行う。
元二十四式の零は次の八式の零と正確に同値である。

```
A1=s*w−f*δ,
A2=s*h−u*δ,
A3=u*w−f*h,
A4=f*P+s*L,
A5=P*w+δ*L,
A6=u*P+s*r,
A7=h*L−r*w,
A8=−(P*h+δ*r)+2*(u*L−f*r).
```

元式への正確な復元は

```
g2=g3=A1, g6=−6A1,
g10=2A4+A2, g11=−2A4−A2,
g12=−2A5−2A3, g14=A5+A3, g15=A5+3A3,
g16=A6, g17=A8, g18=A7,
g21=6A1, g23=−4A1, g24=A2,
```

残g1,g4,g5,g7,g8,g9,g13,g19,g20,g22は零。
逆に A1=g2,A2=g24,A4=(g10−g24)/2、
A3=(g15−g14)/2、A5=(3g14−g15)/2、A6=g16,A7=g18,A8=g17。
従って source parameter を一つも反転しない定数係数の両方向である。
旧 R2=fu はこの正しい系に含まれない。

## 3. 全被約支持は3×3 rank≤1 coneそのもの

次を置く。

```
M = [[s,f,u],
     [δ,w,h],
     [P,−L,−r]].
```

A1,A2,−A3,−A4,−A6,−A5,A7 は M の七つの2minorである。
残二minorを x=uL−fr（rows1,3／cols2,3）、y=−Ph−δr（rows2,3／cols1,3）とすると、
A8=y+2x。

第一row非零なら、上二rowの全minor零から第二row=k*(第一row)となる（k=0も含む）。
既知 rows1,3／cols1,3 minorが零なので y=k*0=0、A8からx=0。
第一row零ならx=0で、A8からy=0。
したがって全九2minorが零、rank M≤1である。
この二場合は全 pointを尽くし、s、δ、P等のgeneric開を仮定しない。
逆にrank≤1なら全minorもA8も零。
従って **全 reduced supportは正確にrank≤1 M**。
非被約scheme idealのradical性や両 missing minorのliteral ideal membershipは主張しない。

## 4. 全 parameterのactual map

L6は全六parameter `(P0,A0,S0,Q0,V0,R0)` に対し

`r=(P0,A0,S0,Q0,V0,0,R0,0,0,−V0,A0,S0)`。

Segre五は全六parameter `(A,B,C,D,E,F)` の外積

`M=(A,B,C)^t*(D,E,F)`

から、元順序で

```
r=(C*D+2*A*F, 2*A*E, A*D, −C*E−B*F, 0, A*E,
   −C*F, B*F, A*F, B*E, 2*A*E+B*D, 2*A*D).
```

任意の rank≤1 matrixは C上この外積で表せ、rank零も含む。
C*の一冗長性があるので coneはdimension五だが、dimensionをcoverの代用には使わない。
元 sourceへの polynomial map は全 parameter上で定義され、inverseは上記rank-one因数分解の存在で尽くされる。
旧反例は `(A,B,C,D,E,F)=(1,0,−Q,0,1,1)` へ入り、fu非零を実際に保持する。

L6とこの全Segre五を合わせ、元十二Poisson sourceの全被約支持を覆う。
旧F/U/S四coneはこのSegre五内のsubfamiliesであり、全源を覆うとはしない。

次はこの **二whole源**だけについて元actualtensor／全350／全126Pfを有限に確認する。
L6の既正しい個別係数結果は再利用可能だが、新Segre五の因子を旧subfamilyから推論しない。
新結果と別作者reviewが揃うまで、whole被約既約quintic空性は未証明である。
