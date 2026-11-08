# K211/K31の全Ore module：simple normal rootはactualDの因子

2026-11-08、PRIMARY G、producer純証明・独立review前。
新CAS、source probe、matrix、GB又はjetは無い。
全K211/K31 strict Qと全quadratic Wに対するactual factorを証明する。
K4/K0へsimple-root論法を外挿しない。

## 1. 全strict symplectic Qに共通するmodule恒等式

full sourceをP=Q(y)+∂c∧W(y)、y=(a,b,u,v)とする。
Qはstrict quadratic Poisson、K4=PfQは非零normal binary quartic。
全projective JacobiはE_y∧L_WQ=0。

normal chartv1ではpi=projectivized Q、Z_i=Q^{i,v}|v1として

`P_chart=pi+∂c∧(Y+cZ)`, Y_i=W_i−y_iW_v。

Qをradial logarithmic coordinateとprojective chartへ分解すると、strict Jacobiから
L_Zpi0で、piはrank二、Zはimage piへ入らない。
inverse symplectic formω_Qはhomogeneous weight零。
Cartan formulaからi_Eω_Qはclosedであり、radial quotientへdescendする。
符号をβ(Z)=1へ固定すると、piを殺すclosed rational one-form βが得られる。
これがradial quotientの正確なcosymplectic pairである。

local primitive r of βとZ-flowでZ=∂rを取り、残二coordsをx,zとすれば
pi=ρ(x,z)∂x∧∂z、ρはr非依存。
元module L_Ypi+Y∧Z0の(r,x),(r,z)成分から

`Y=H Z+pi#dH`, H=β(Y)                               (1)

が従い、残成分はHamiltonian恒等式になる。
逆も成立する。
この一時local analytic座標をactual PGLと呼ぶのではなく、
rational identity(1)を元座標へ戻して使う。
全Wが含まれ、YからWへのkernelはEuler redundancyだけ。

actualDは

`D=K4*(c+Hhom)`。                                      (2)

Hhomはdegree一で、K4 Hhomはpolynomial F5。
normal linear factorLがK4へ一回だけ現れれば、HhomのL-pole orderは高々一。
chartでのresidue h(a,b)はpolynomialである。
K4がnormal variablesだけのpolynomialで、他normal factorsはgeneric chartpointでunitだからである。
residueにa/bの未許可分母を導入しない。

## 2. 全K211：L=u−vをactual因子にする

whole canonical sourceは全t∈Cと全C2に対し

```
U=(u²,uv),
V=(u²+(t−2)uv,−2uv+(t+1)v²),
A1=−u+(2−t)v, B1=u.
```

K4はnonzero scalar倍のu²v(u−v)。
chart v1でliteral projective coefficientsは

```
pi_au=0,
pi_bu=3u(u−1),
pi_ab=3a(1−u)+C2(u,1),
Z=(u,−2u+t+1,0).
```

simple root u1でH=h(a,b)/(u−1)+regularと書く。
(1)のa成分は

`Y_a=H*u+pi_ab H_b`。

pole coefficientは

`h+C2(1,1)h_b`。

Yはglobal quadratic Wから来るpolynomial chart vector fieldなので、これが零でなければならない。
nonzero constant kとconstant dに対し、k+d∂bは有限polynomials上可逆：
degree in bのhighest coefficientから順にzeroとなるか、finite nilpotent geometric seriesでinverseを書ける。
ここではk1だからh0。

従ってHhomはL=u−vにpoleを持たず、(2)からactualDはLで割れる。
t、C2、Wの全parameterを保持したfactorであり、irreducible quinticを持つ全K211源は空。

## 3. 全K31：simple rootvをactual因子にする

whole sourceはr²+st0の全nilpotent matrix、C零も含め

```
U=(u²+r uv,s v²),
V=(t uv,uv−r v²),
A1=(r−2t)v, B1=u+(2r−s)v.
```

K4はnonzero scalar倍のu³v。
chart u1でnormal ratioをvとし、literal coefficientsは

```
pi_av=−v+(s−r)v²,
pi_bv=v−(r+t)v²,
pi_ab=(r−t)a v+(r−s)b v+C2(1,v),
Z=(1+r v,t v,0).
```

simple rootv0でH=h(a,b)/v+regular。
(1)のa成分のpole coefficientは

`h+C2(1,0)h_b+h=2h+C2(1,0)h_b`。

最後のhはpi_av H_v=(-v)*(-h/v²)からである。
2+d∂bも全有限polynomials上可逆だからh0。
Hhomはvにpoleを持たず、Dはvで割れる。
projective nilconeの全P1、C零、全magnetic C2及び全Wを同時に覆う。
従って全K31源もirreducible quinticを持たない。

## 4. 正確な結論と残余

先行strict Ore bridgeと全binary-Q classificationの適用hypothesesの下で、
K211/K31は**全actual factorで空**。
K22は別whole hand module proofにより一actual projective orbit、G像≤24。
残るnonzero quartic branchはK4a/K4b、及びK零の非constant-factor源。
今回のsimple-pole coefficientをmultiple poleへ流用しない。
全concurrent G、full E recognition、ambient component数にはまだ未解決部分がある。
