# Jordan (2,1)、(L,M)=(h,x) の全六源を半単純族へ覆う実変形

2026-11-08、PRIMARY `E_reduced_nonnormal`。新成分を作らず、既存 MPS40 の閉和への全源被覆を示す候補である。
以下の手導出を、元の全六源・全150係数・全350式と照合する一つの有限 job を実行する。
独立 review は未完。以前の指定点の接次元、法線、局所 jet を用いない。

## 1. 元の全六係数と全逆

`TWO_HYPERPLANES_CUBIC_SUPPORT.md` §6 および
`code/data/classification_two_hyperplanes_cubic_nilpotent_sources.json` の
`j21_h_x_irreducible` は、全45 vertical 係数から全 Jacobi と三次既約性を保って得た必要十分な全六源である。
自由変数は `n0,n5,n7,n8,n10,n12`、既約性には `n5!=0` が必要である。
斉次座標は `(h,x,u,v,w)`、h 行は零、cofactor は

```
V0=v∂u+h∂v+x∂w.
```

`A=n12, B=n12*n8+n7, C=n0+n10*n12` と置き、次の三角 PGL を取る。

```
U=u-Aw-Bx-Ch,  V=v-Ax,  W=w-n10*h-n8*x.
```

det は 1。逆は `w=W+n10*h+n8*x`, `v=V+Ax`, `u=U+Aw+Bx+Ch` であり、
source parameter を一つも反転しない。この変更は V0 を exact に保つ。
`s=n5, k=n8` とすると変更後の全括弧は正確に

```
{x,u}=xv,              {x,v}=hx,            {x,w}=x²,
{u,v}=v²/2-hu,         {u,w}=s h²+k xv-vw,  {v,w}=k hx-hw.       (1)
```

実退化五次式は

```
D0=hx[-s h³+x(v²/2-hu)].                                     (2)
```

s が非零なら三次式は u について一次で、係数 `-hx` と残部
`xv²/2-s h³` は互いに素なので既約である。s=0 は h または x の追加因子を持つ補集合で、
元の既約三次枝ではないが以下の同じ正則曲線はこれも保持する。
全源の四つの余分な自由係数を、指定点へ置いて落とした議論ではない。

参考として `x->x/s,w->w/s` は (1) を s=1 へ送り、
`t^-1 diag(1,t^-2,t²,t,t^-1)_*` は s を保って k を tk へ送る。
従って k 非零は一つの軌道、k=0 はその実多項式境界である。
しかし下記曲線は全 s,k を直接覆うため、この軌道論や source 次元を被覆の仮定にはしない。

## 2. 全三 parameter の真の半単純親

親座標 `(H,X,U,V,W)` に対して

```
{X,U}=3XU,  {X,V}=XV,  {X,W}=-XW,
{U,V}=aUV,  {U,W}=-bV²-aUW,  {V,W}=A H²+aVW.                 (3)
```

全自由 parameter `(A,b,a)` に対して strict Poisson、実五次式は

```
H X[3A H²U+bV³+3aUVW].                                     (4)
```

これは既受理 `TWO_HYPERPLANES_CUBIC_SUPPORT.md` §3.2 の
半単純固有値 `(3,1,-1)`, j=0 の全 affine 源である。
`A=b=0,a=1` は `3HXUVW` の rank 四 SNC 点を与える。
全 A3 源の像閉包は既約なので、それを含む一つの ambient 成分はこの SNC 点を含む。
既受理の全 SNC-meeting component 尽くし定理により、その ambient 成分は MPS40 の一つである。
従って (3) の全源が既知 MPS40 の有限閉和に入る。
一点の既知帰属を無断で別の族へ伝播していない。

## 3. 全 pole を消す actual source 曲線

`e!=0` で、真の親の parameter と座標を次の通り代入する。

```
H=h, X=x,
U=w+x/(3e),
V=u+v/(2e)+h/(2e²),
W=v-h/e,
a=1-3ke,  A=a/(2e³),  b=-8s e⁴,
tensor scale=e.                                             (5)
```

座標変更の det は 1、逆は

```
w=U-X/(3e),  v=W+H/e,  u=V-W/(2e)-H/e².
```

e のみを反転し、s,k を反転しない。積 `A=a/(2e³)` を親 parameter の真の代入のまま保持する。
特に `Q=v²/2-hu` とすれば

```
V*W=uv+Q/e-h²/(2e³),
A*h²+a*V*W=a*(uv+Q/e).
```

これは最も大きい pole の actual 相殺である。全 tensor の正則表示は

```
{x,u}=xv+e*xu,        {x,v}=hx-e*xv,       {x,w}=x²+3e*xw,
{u,v}=a*(Q+e*uv),
{v,w}=e*b*V²+e*a*wv-a*hw+(a-1)*xv/3-(a-1)*hx/(3e),
{u,w}=-a*e*uw-a*vw+(1-a)*xu/3+(1-a)*xv/(3e)-b*V²/2.        (6)
```

最後二式に (5) の a,b,V を代入すると全係数が QQ[e,s,k] の多項式である。
e=0 の六括弧は厳密に (1) になる。全 k=0、s=0、他の全係数零の境界を同じ式に保持する。
各 e 非零は (3) の actual source 点への可逆 PGL と非零 tensor scale であり、
e=0 は pole-free な実 Poisson tensor の端点だから MPS40 の閉和に入る。

## 4. 全六源の結論と正確な補集合

有限 job で元の全 source binding・三角 normalizer/逆・全六括弧・全150係数・
親/実曲線の全350および strict10 Jacobi・実 Pfaffian と (4) の実 pullback を照合する。
成功し独立 review が受理されれば、元の全六 source の全 parameter と PGL 飽和閉包が MPS40 に入る。
四つの余分な係数の消去は全逆のある三角変更なので、generic parameter だけの主張ではない。
新 component は零、ambient smoothness/local uniqueness は不要である。

`(L,M)=(h+x,h)` の p=0 全五源、`(0,h+x)` の u 依存全六源および n6=0 補集合、
他の nilpotent 源と平面 1+1+3 の残三配置はこの結果の対象でなく、未解決のまま残す。
