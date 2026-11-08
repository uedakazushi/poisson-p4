# contact31_chordの全六磁気源：二conicのSNC族からの実衝突

2026-11-07、PRIMARY `E_reduced_nonnormal`。計算前に全源と実変形を固定する。
成分の新規追加は目的とせず、既知MPS有限閉和への全源帰属を狙う。

## 1. SNC点を含む既約な実Poisson源

親座標を(H,A,B,C,D)とし、p,q,k,u,vを任意の複素parameterとする。
N1=A∂A−B∂B、N2=C∂C−D∂Dと置き、

```
P=p(AB+uH²)∂A∧∂B+q(CD+vH²)∂C∧∂D+k N1∧N2
```

を使う。二つの二次元Poisson tensorは異なる変数に依存して可換であり、
N1,N2は互いに可換で、それぞれ両tensorを保存する。従って全五parameterでstrict Poissonである。
実affine Pfaffianはpq(AB+uH²)(CD+vH²)、射影因子はさらにH倍となる。
u=v=0,pq≠0は五つの独立超平面H,A,B,C,DのSNC点である。

親parameter空間A5は既約なので、全像閉包を含む全Poisson成分を一つ選ぶと、
その成分はこのSNC点を含む。既受理のMPS Theorem1.4の
「SNC log-symplectic locusと交わる全成分は既知40個」というscopeにより、
親の全像はMPS40の有限閉和へ入る。特定のT番号又はSNC点の成分一意性は不要である。

## 2. 一つの実PGL衝突と全periodの積関係

target座標を(h,x,y,z,w)とする。u=v=−1とし、

```
(H,A,B,C,D)=(h,x,y,x+epsilon z,y+h+epsilon w),
p=−1+epsilon s, q=1, k=−1+epsilon t
```

を親へ代入してtensorをepsilon倍する。座標行列式はepsilon²である。
epsilon≠0では実PGL変換と非零scalarなので、全点が上の既知閉和に入る。
parentのp,q,k及びedge u,vは元の実源の関係のままであり、
独立Taylor product係数へ緩めたものではない。

Q1=xy−h²と置くと、正則極限の括弧は

```
{x,y}=0,
{x,z}=−x², {y,z}=h²,
{x,w}=2xy+hx−h², {y,w}=−y²−hy,
{z,w}=(2y+h)z+s Q1+t hx.
```

分子のepsilon^(-1)項は、
`(p+q)Q1+(q+k)hx=epsilon(sQ1+t hx)` により全て消える。
全epsilon族は斉次二次の多項式tensorである。

## 3. 任意の六磁気係数を端点に戻す

accepted contact31_chord源は
`reports/CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_20261107.json` にあり、
括弧は上記U=(−x²,h²)、Vの−5/7倍、fiber一次項も−5/7倍で、磁気項は

```
M=m0h²+m1hx+m2hy+m3x²+m4xy+m5y²
```

である。c=−5/7とする。
極限前の実族に、epsilonに依存しないfiber translation
`z'=z+ell, w'=w+kappa` と `W'=c w'` を施す。
次の全六parameterの式を固定する。

```
s=−7m4/5, t=−7m1/5,
a=7m2/10−7m5/15, b=7(m0+m4)/5−a,
d=7m5/15, e=−7m3/5,
ell=a h+b x+d y, kappa=e x.
```

U,Vとfiber構造係数から、translationの磁気変化は正確に

```
−(2y+h)ell−U(kappa)+V(ell)
=−(a+b)h²−2(a+d)hy−3d y²+e x².
```

これをsQ1+t hxへ加えてc倍すると、全係数がMに一致する。
baseは不変なので、translationの逆はell,kappaを引くだけであり、
最後のfiber scaleの逆も定数−7/5である。
磁気係数を一つも反転せず、全零磁気項を含む全境界を直接扱う。

したがって変換後の実族Pi_epsilon(m)は、epsilon=0でaccepted源の全六括弧に一致する。
この全係数一致、全strict10/元350Jacobi、Pfaffian、正則性を一つの有限jobで確認する。
陽性なら閉性により任意のmで端点も既知閉和に入り、全PGL飽和が被覆される。
局所一意性、接上界、任意の高次jetは使わない。

## 4. 保持する補集合

唯一の非定数反転epsilonは親をPGL変換として解釈するためだけで、その零fiberが対象である。
pが途中で零となるrank低下点も、親の全既約源内の点として保持する。
全mを代入できる多項式族なので、source parameter分母の補集合は無い。
もう一つの残平面源contact31_tangent_at_simple_intersectionへの一般化は主張しない。
