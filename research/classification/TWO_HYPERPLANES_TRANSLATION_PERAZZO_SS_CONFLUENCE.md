# Perazzo三次・二超平面の全源：半単純(2,1,−1)族からの実変形

2026-11-08、PRIMARY `E_reduced_nonnormal`。全source有限化後のactual parent橋。
七法線のjetや別canonical点を選ばず、同じkappa=1モデルを既知閉和へ入れる。

## 1. SNCを含む全三parameterの真の親

座標(H,X,U,V,W)、自由parameter A,L,aに対し、H行を零、残りを

```
{X,U}=2XU, {X,V}=XV, {X,W}=−XW,
{U,V}=aUV,
{U,W}=−LH V−(a+2)UW,
{V,W}=AHX+aVW
```

とする。これは `TWO_HYPERPLANES_CUBIC_SUPPORT.md` §3.2の
semisimple(2,1,−1)、j=1,ell=0の全sourceそのもので、全A,L,aでstrict Poissonである。
実五次は

```
H X[2AHXU+LH V²+2(a+1)UVW].
```

A=L=0,a=1は五つの独立超平面のSNC点。
親A3は既約なので、その全像を含むambient成分はこのSNC点に会い、
採用済みの全SNC成分尽くし定理からMPS40の一つとなる。
従って親の全像は既知MPS有限閉和内にある。
単なるSNC一点の所属から曲線の所属を逆推論しない。

## 2. 導出されたLaurent profileと全pole相殺

internal座標を(h,x,u,v,w)とし、epsilon≠0で真のaffine/PGL shear

```
H=h, X=x,
U=u+h/epsilon,
V=v+x/epsilon,
W=w+(h−x)/epsilon
```

を行う。detは1、逆も同じ一階poleだけを持つ。
tensorをepsilon倍し、親の**真のparameter**を

```
a=−2+2kappa*epsilon,
A=epsilon^(−2)−2kappa*epsilon^(−1),
L=−2epsilon^(−1)+4kappa
```

へ代入する。a+2とA,Lの積を独立Taylor係数へ緩めていない。
このprofileはshiftの三linear forms h,x,h−xと、元の全括弧の負冪相殺から直接得た。
inverse shearの一階poleと二次係数から、上の表示に現れるpoleを全て文字通り検査できる。
arcの次数を順に増やす探索ではない。

最後にepsilonに依存しない可逆PGL

```
G:(h,x,u,v,w)→(h,−2x,−2v,u/2,u/2−v+w), det G=−2
```

を施し、最終座標を再び(h,x,u,v,w)と書く。
得るPi_epsilon(kappa)は全括弧でepsilonの負冪が消えた、斉次二次の多項式tensorとなる。

## 3. 同じcanonicalモデルへの全六端点恒等式

epsilon=0では全kappaについて

```
{x,u}=x², {x,v}=hx, {x,w}=0,
{u,v}=hu−kappa hx,
{u,w}=−hu+xv−xw+kappa x(x−h),
{v,w}=−hu+xv+kappa h(x−h).
```

これは `TWO_HYPERPLANES_TRANSLATION_NONPLANE_CANONICAL_MODELS.md` §2と
completed nonplane reportの同じPi_kappaである。
kappa=1で全D96/orbit47を計算した元の点そのものに一致する。
新たな点のrankを調べていない。

F0=(h−x)(hu−xv)+hxwと置くと、真の全epsilon五次は厳密に

```
D_epsilon=h*x*(1−2kappa*epsilon)*[
 F0+epsilon*(hu(v+w−u)+2xv(w−v))
    +epsilon²*uv(2w−2v−u)].
```

従って端点は全kappaでhxF0≠0を保つ。
このDだけの一致で帰属を推測せず、全150tensor係数、全strict10/元350 Jacobi、
全Pfaffian、shear/Gの逆、全負冪相殺を一つの有限jobで確認する。
kappaは一度も反転せず、零kappaも同じ正則端点として含む。

## 4. 全七parameter源への最弱の伝播

epsilon≠0は親の真の点に可逆PGLと非零scalarを施したものなのでMPS有限閉和へ入る。
全括弧がepsilon=0で正則であるため、任意のkappaの端点も同じ閉和内にある。

`TWO_HYPERPLANES_TRANSLATION_NONPLANE_ORBIT_CLOSURES.md` の
13次元flag base＋9次元固定N fibre上界、実22orbit下界を使えば、
非平面L=h−xの全七parameter源のPGL飽和閉包はkappa=1 orbit閉包に等しい。
したがってこのcanonical点の帰属が、その全source、全normalizer、全parameter境界へ伝播する。
局所一意性、七法線の支持、formal IFTの高次は不要である。
この全source orbit橋及び新arcの独立review後にのみ、当該E枝を受理済み被覆へ移す。

L=hのessential四変数源、平面A=B=0、他の1+1+3源又は全Eは別に残す。
