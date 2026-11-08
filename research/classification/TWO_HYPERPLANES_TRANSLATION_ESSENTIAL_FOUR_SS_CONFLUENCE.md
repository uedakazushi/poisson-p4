# essential四三次・二超平面の全源：半単純(3,1,−1)族からの実変形

2026-11-08、PRIMARY `E_reduced_nonnormal`。Perazzo源と混ぜず、L=hの全源を扱う別のactual親。
七法線、高次jet、別kappa代表を用いない。

## 真の全三parameter親とSNC

座標(H,X,U,V,W)、自由parameter A,b,aに対して

```
{X,U}=3XU, {X,V}=XV, {X,W}=−XW,
{U,V}=aUV,
{U,W}=−b V²−(a+3)UW,
{V,W}=AHX+aVW
```

とする。H行は零である。
これは `TWO_HYPERPLANES_CUBIC_SUPPORT.md` §3.2の
semisimple(3,1,−1)、j=1の全源で、全A,b,aについてstrict Poissonである。
実五次は

```
H X[3AHXU+b V³+3(a+1)UVW].
```

A=b=0,a=1は実因子6HXUVWのSNC点。
親A3は既約なので、全像を含む一つのambient成分はこのSNC点にも会い、
全SNC成分尽くし定理によりMPS40の一つである。
従って親の全像が既知MPS有限閉和に入る。

## 全poleを消す一つの実source曲線

internal座標(h,x,u,v,w)上で

```
H=h, X=x,
U=u+h/epsilon,
V=v+x/epsilon,
W=w+h/epsilon,
a=−3+6kappa*epsilon,
A=2epsilon^(−2)−6kappa*epsilon^(−1),
b=−6epsilon
```

とし、tensorをepsilon倍する。shearのdetは1、逆も一階poleだけを持つ。
a+3及びA,bは全て元の実parameterのまま代入し、係数の積を緩和しない。
最後に定数PGL

```
G:(h,x,u,v,w)→(h,x,v,−w,u/6+w/2), det G=−1/6
```

を施す。
全epsilon負冪が六括弧から消え、全kappaに定義された斉次二次の多項式tensorを得る。

端点は

```
{x,u}=x², {x,v}=hx, {x,w}=0,
{u,v}=3hu−2xv,
{u,w}=−hu+xv−kappa hx,
{v,w}=x²−kappa h².
```

これは同じcanonical model note§1及びcompleted source reportのPi_kappaである。
元の全八parameter源の固定点からの既知変更v_new=v_old−2wを使っただけで、
新しいPoisson点やrank評価を選んでいない。
その逆変更を端点kappa=1へ戻す全150一致も検証する。
実端点因子は全kappaで

```
hx(h²u−hxv+x³) != 0.
```

## 全源へ必要な最弱橋

親の全Jacobi/SNC/Pfaffian、shearとGの全逆、負冪相殺、端点全150、
元350及びstrict10 Jacobi、全epsilonのactual Pfaffianを一つの有限jobで確認する。
epsilon≠0は既知親の真の点、epsilon=0は正則端点だから全kappaのcanonical点が既知閉和へ入る。
kappaを一度も反転せず、kappa=0も同じ実族で含む。

13次元cofactor/二Hflag base＋9次元exactN fibre、実22orbit下界の独立reviewが受理されれば、
`TWO_HYPERPLANES_TRANSLATION_NONPLANE_ORBIT_CLOSURES.md` により
L=h全八parameter源のPGL飽和閉包はこのkappa=1 orbit閉包そのものである。
よって全parameter、n17非零の全域、H/x交換したchart、その他のnormalizer、全source境界へ被覆が伝播する。

平面A=B=0、他のnilpotent/cofactor型及び全E被覆は別に残す。
