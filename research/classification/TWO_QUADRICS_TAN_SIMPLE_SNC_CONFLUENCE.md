# contact31_tangent_at_simple_intersectionの全磁気源：coupled SNC族の実衝突

2026-11-07、PRIMARY `E_reduced_nonnormal`。有限検証前に実族と全六磁気の逆写像を固定する。
他の固定profile・別anchor又は高次jetを使わない。

## 1. 積関係を保つ全五parameter親

親座標は(H,A,B,C,D)。p,q,a,u,vを自由とし、b=p−a、
F=pAB+uH²、G=qCD+vFと置く。

```
P=F∂A∧∂B+G∂C∧∂D
  +(A∂A−B∂B)∧[a C∂C+b D∂D].
```

これは全五parameterでstrict Poissonである。
例えば(A,C,D)のJacobiは
`vpAF+pqACD−(a+b)A(qCD+vF)=0`、他の三式も同じ関係で零となる。
b=p−aとG=qCD+vFを元のまま保つ必要があり、vFの各係数を独立化していない。
実affine PfaffianはFGで、u=v=0では射影因子pqHABCDとなる。
従って全親A5は既約でSNC点を含む。
全像を含むambient成分はそのSNC点を含むので、既受理MPS SNC成分尽くし定理により
全親はMPS40の有限閉和へ入る。SNC点が既知閉和に属するだけという推論は使わない。

## 2. 真のPGL衝突と全pole消去

targetを(h,x,y,z,w)とし、

```
(H,A,B,C,D)=(h,y+h,x−y+h,x+epsilon z,−y+epsilon w),
p=1, u=−1, v=−1/3,
q=−1/3+epsilon q1, a=2/3+epsilon a1, b=1−a
```

を代入してtensorをepsilon倍する。座標行列式は−epsilon²。
epsilon≠0は実PGLと非零scalarなので、既知MPS閉和内にある。
二normal linear formsはC−A−B+2H、D+A−Hである。
その括弧をplaneへ制限したepsilon定数項は零であり、epsilon一次項は正確に

```
−q1 xy−a1(hx+2y²).
```

従ってtensor epsilon後の唯一あり得るepsilon^(-1) poleは完全に消える。
端点anchorは

```
U=[2x(2y−x)/3, −(hx+xy−3y²)/3],
V=[(y²−2xy−3hx)/3, −y(y+h)/3].
```

そのcommutatorはαU+βV、α=(2h+y)/3、β=2y/3で、fiber括弧は

```
{z,w}=−αz−βw−q1xy−a1(hx+2y²).
```

全実族の正則性、四anchor、fiber一次項及び全Jacobiは一つの有限全係数jobで確認する。

## 3. 全六磁気sourceへの可逆移送

accepted sourceは
`CLASSIFICATION_TWO_QUADRICS_PLANE_SOURCES_20261107.json` の
`contact31_tangent_at_simple_intersection`。
baseの二conicはQ1=hx−y²、Q2=hx+xy−y²である。
定数fiber変更

```
(Z,W)^T=T(z,w)^T,
T=[[-5/4,0],[25/24,5/2]], det T=−25/8
```

により上のU,Vはaccepted二anchorへ正確に移る。
目標磁気項M=m0h²+m1hx+m2hy+m3x²+m4xy+m5y²に対し、
まずfiber scale前の係数ni=−8mi/25を置く。
次の六parameterは全miのlinear formである。

```
a1=−(n1+n5)/3,
q1=−n4−n3−n5/3+2n1/3,
A=3n0/2, B=n5−2n1, C=3n2−A, E=3n3/2.
```

epsilonに依存しないfiber translationを
`z'=z+ell, w'=w+kappa`、ell=Ah+Bx+Cy、kappa=Exとし、その後にTを施す。
translationの磁気変化は厳密に

```
αell+βkappa−U(kappa)+V(ell)
=2Ah²/3+(A+C)hy/3−Bhx/3
 +(−B−2E)xy/3+By²/3+2Ex²/3.
```

これを−q1xy−a1(hx+2y²)へ加えるとΣni monomial_iとなり、
最後のdet T倍により全磁気項がMと一致する。
translationの像はspan{h²,hy,−hx−xy+y²,x²−xy}でrank4、
親の残二period −xy,−hx−2y²はその商でrank2である。
したがって全六磁気は同時に独立であり、各零係数境界も残る。

translationの逆はell,kappaを引き、fiberの逆はT^(-1)である。
六個のmから(q1,a1,A,B,C,E)への行列は非零定数detを持つ。
source parameterの反転、period-edge productの緩和、磁気一点への固定をしない。

## 4. 最弱の結論

上の真の全epsilon族が全150端点、全350/strict10、実Pfaffianで確認されれば、
一般epsilonは既知MPS有限閉和、正則端点はaccepted全六磁気源である。
閉性によりこの源の全点と全PGL飽和が被覆される。
局所一意性、新成分、全schemeの滑らかさ、三次以上の支持は要求しない。
唯一の非定数反転epsilonは親PGLの解釈に用い、その零fiberそのものが対象である。
