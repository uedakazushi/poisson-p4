# Concurrent位数三・K22 Ore枝：全源は一actual軌道、G上界24

2026-11-08、PRIMARY G、producer純証明・独立review前。
新CAS、source probe、module matrix、GB又はjetは無い。
先行全row-pair還元の二fixed Qを全Wのまま解く。
結論は全K22枝のprojective PGL5像次元≤24であり、ambient componenthood又は既成分E所属は主張しない。

## 1. 先行全source還元を適用する正確な対象

actualDは被約既約quintic、rank二planeのzero cubic fは三reduced concurrent lines、normal位数三。
全Ore liftのquartic coefficient K4=Pf(Q)がpartition22とする。
`PLANE_CONCURRENT_M3_K22_REDUCED_CUBIC_REDUCTION.md`
SHA5bd2a0efe8833cdcd6dc5f020cf6ff948e2c76dbacbc3e8cd55c179fe358f2f0
の全p/q zero complements・actual factorを使う。
残源は二fixed γ=0,1について

```
Qab=−a v+b u+γuv,
Qav=Qbu=uv, Qau=Qbv=Quv=0,
P=Q+∂c∧W(a,b,u,v), E_y∧L_WQ=0.
```

full quadratic Wの全first/higher coefficientsを保持する。
Wはnormal ideal(u,v)を保ち、base quadratic partはEuler gauge後(f_b,−f_a)である。
P0 first symbolだけのsliceから全Wを推論していない。

## 2. 全moduleは一つの有理Hで一意に決まる

chart v1で

```
pi_ab=−a+γu, pi_au=−u², pi_bu=u,
Z=u∂a,
P_chart=pi+∂c∧(Y+cZ),
Y_i=W_i−y_i W_v, y_i=(a,b,u).
```

元Jacobiは `L_Ypi+Y∧Z=0`。
closed rational one-form

`β=d(a/u+b)+γ du/u`

はpiを殺し、β(Z)=1。
local analytic r=a/u+b+γloguを用いれば、pi=u∂b∧∂u、Z=∂rとなる。
この一時的なanalytic座標は計算のためだけで、ambient PGL変換とは呼ばない。

H=β(Y)と置く。
moduleの(r,b),(r,u)成分はY_b=uH_u,Y_u=−uH_bを強制し、
(b,u)成分はその恒等式になる。
従って元のrational coordinatesで全sourceは正確に

```
Y=H Z+pi#dH,
Ya=uH+(−a+γu)H_b−u²H_u,
Yb=(a−γu)H_a+uH_u,
Yu=u²H_a−uH_b.                                       (1)
```

逆に任意Hに対しこれらはmoduleを満たす。
YからWへのkernelはW=ell(y)E_yであり、
Pの差はactual Euler redundancyなので同じprojective tensor。
**Hが決まれば全ambient tensor classも一意**で、未読higher extensionの自由parameterは無い。

## 3. ActualDと全HのLaurent次数上界

K4=Pf(Q)=u²v²。
actualDはoverall scalarを固定すれば

`D=u²v²*(c+Hhom(a,b,u,v))`。                          (2)

Hhomはdegree一rational homogeneous functionで、chart v1での値がHである。
これはE∧P²/2又はpi∧(Y+cZ)のliteral計算で従う。
F5=D−c u²v²はpolynomialでnormal order≥三。
従ってHhomのdenominatorは高々u²v²、numeratorは全degree五/normal order≥三。

chartでは

`H=Σ_(k=−2)^3 h_k(a,b)u^k`,

deg h_k≤二、さらにdeg h_2≤一、h_3はconstant。
全二十八monomialを初めから許す。
無限jet又はformal power seriesを用いない。

## 4. u0及びv0の全poleを比較する

(1)がpolynomial Yを与える必要条件を全negative u coefficientsへ課す。
Yuのu^-1から(h_-2)_b=0。
Ybのu^-2からa(h_-2)_a−2h_-2=0。
よってh_-2=A a²。

Yaのu^-1からa(h_-1)_b=3A a²、従ってh_-1=3Aab+g(a)。
Ybのu^-1から

`a g'(a)−g(a)=2γA a`。

gはdegree≤二のpolynomialなので、全solutionは

`γA=0, g=B a`。

同じactual swap(a,b,u,v)→(b,a,v,u)ではγ→−γとなる。
Hhomのdegree一twistを保持してchart u1の全negative v coefficientsへ同じ比較を適用する。
これにより全remaining h_kが尽くされ、必要な全Hは正確に

```
Hhom = A*(a²v/u²+3ab/u)
      +A'*(b²u/v²+3ab/v)
      +B*a v/u+B'*b u/v+ell(a,b,u,v),                  (3)
γA=γA'=0,
```

ellは全linear formである。
例えばv^-2 coefficientはA'b²、v^-1 coefficientは3A'ab+B'bでなければならず、
その他のdegree二/一/零項を全て比較した。
u0だけを調べてv0やfull P1/P2を捨てていない。
全sourceが(3)へ入るという必要方向だけでも以下の上界には十分であり、
未確認のformal extensionの一意性は仮定していない。

## 5. γ一は位数三源を持たない

γ1ならA=A'0。
(2)(3)のF5はnormal order四以上であり、G3零。
今回のnormal order三に反するためこの二fixed moduleの一方は**全空**。
B/B'/ellを零へ選んだpointの話ではない。

## 6. γ零の全reduced sourceと真の正常化

γ0とする。
(1)のhighest base degree部分を、唯一のbase-gradient gaugeへ戻すと

`f=2ab(A b−A'a)`。                                    (4)

具体的にはA/A'部分のpolynomial W0は

```
Wa=−2A'a²+4Aab,
Wb=4A'ab−2Ab²,
Wu=−2A'a u+A b u−A a v,
Wv=−A'b u+A'a v−2A b v.
```

Wa=f_b,Wb=−f_aである。
従ってf reducedは正確にAA'≠0。
AA'0の補集合を無根拠に捨てず、stated reduced-R hypothesisが除くことを明示する。

ellはactual c'=c+ell(y)によって(2)のHhomから除ける。
Qは変わらず、Hからのmodule一意性によりfull tensorも同じsourceへ戻る。

actual a'=a+αu,b'=b+δvはγ0のQを正確に保つ。
(3)へ代入するとlinear ellの変化を除き

```
B →B−A(2α+3δ),
B'→B'−A'(3α+2δ).
```

determinantは−5AA'≠0なので、全B,B'を一意に零へ戻せる。
zero/nonzero higher directionsをchosen sliceへ制限せず、fullsourceのactual invertible shearで除いた。
inverseはminus shearである。

最後にpaired normal/base torus

`(a,u)→s(a,u), (b,v)→t(b,v), c→c`

はQを保ち、(A,A')を(A/t,A'/s)へ動かす。
t=A,s=A'で両方一へ戻る。
全scalesはnonzero、inverseもactual diagonalである。

## 7. Canonical tensorと全G上界

従って全対象は一つのactual projective PGL5 orbitへ入る。
canonical Qはγ0、canonical Wは

```
Wa=−2a²+4ab, Wb=4ab−2b²,
Wu=−2a u+b u−a v,
Wv=−b u+a v−2b v.
```

全tensorはP=Q+∂c∧W。
全actualDは

`D=c u²v²+a²v³+3ab u v²+3ab u²v+b²u³`。

上のrational module恒等式からprojective Jacobiがdense openで零、従ってglobal polynomialにも恒等零。
Dはc一次、coefficient u²v²と残項のgcdは一である。
u0ではa²v³、v0ではb²u³なので、Dはprimitive linear polynomial in cとして既約。
G3非零、canonical f=2ab(b−a)は三distinct concurrent linesである。

PGL5のdimensionは24だから、このwhole orbit及びそのclosureのdimensionは高々24。
これで**全K22 Ore枝はG非優越**。
orbitの正確なdimension、ambient componenthood、全既componentとのE認識は今回の主張ではない。
K211/31/4/K0の残whole modulesとconcurrent G全体は未解決のまま保持する。
