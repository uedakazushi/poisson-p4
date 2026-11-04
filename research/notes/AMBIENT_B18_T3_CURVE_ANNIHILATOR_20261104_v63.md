# Top T3をactual CA平面曲線上のmoduleへ還元：96恒等式と一つのcolon

2026-11-04。既構成の全186相対入力を使い、T3に関する新しい正確な還元を行う。
T3の実際の零・非零、colon、socle、全Dのdepth19/20やnilradicalはまだ計算していない。
新しいGB、source/core再生成、点、jetは使わない。

## 1. 全96本のannihilationはmatrix内の列で明示される

R0²内のfull matrixを前packetのまま使う。T3=Jinf/Q[3]は、六cubic pair c_alphaに
actual CAのpositive odd coefficients j1,j2を掛けた十二列

```
t_(alpha,e)=(j_e*tau_A)*c_alpha
```

で生成される。I0_CAの八generatorをf_iとする。すると

```
f_i*t_(alpha,e)=j_e*((f_i*tau_A)*c_alpha).      (1)
```

右の括弧はQ[3]の実際のC0-columnである。全indexは0-basedで、alpha=0..5、
e=0,1、i=0..7に対し

```
top column=2alpha+e,
Q[3] column=1980+9alpha+(i+1),
cofactor=j_e.
```

具体的にc_alpha=(a,b)なら、左辺は `(nu_A*f_i*j_e*b, f_i*j_e*a)`、
右辺も同じ。96恒等式の両行を、DAGの可換積normal formとして検査し、全pointerと
cofactor nodeを保存した。formal coefficientを有限jetで零とした検査ではない。
従って **I0_CA*T3=0** が全次数で成り立つ。

## 2. Actual CAの九座標と平面曲線

前packetの48unit source DAGにexact rational ADを行い、actual CAのzero odd八本と
normalized positive odd二本について10×9 Jacobianを求めた。
source48 Jacobian determinantは−84934656/625、選択九本

```
z1,...,z7,r,s
 = f18,f19,f28,f29,f60,f67,f72,j46,j95
```

のdeterminantは **128/3**。row63の一次部分は零である。
従ってこれらはR0のactual formal/algebraic coordinates。
一つのpointed57unit inverse系を追加し、

```
g_A(r,s)=f63(X(z=0,r,s))
```

を有限DAGで指定した。coordinate telescopingにより

```
I0_CA=(z1,...,z7,g_A),
R0/I0_CA ≅ C[[r,s]]/(g_A)=C_A.                 (2)
```

g_Aが零なら、full odd ideal I_CA=I0_CA+wJ_CAを七つのzとwr,wsの九本で生成でき、
既受理のactual mu(I_CA)=10に反する。よってg_A!=0。
従ってC_Aは一次元Cohen–Macaulay plane curveであり、(1)からT3はその上の有限moduleとなる。
actual g_Aのfactorizationや展開した全係数はここでは求めない。

## 3. Multiplicity3は旧analytic modelの座標同一視なしに従う

同じactual coordinatesで、full conductor quotientの共通G19横断環は正確に

```
Ccore=C[[w,r,s]]/(g_A,wr,ws)
      ≅ C[[w]] fiber_product_C C_A.
```

受理済みの群作用による同時product decompositionにより、conductorのmultiplicity4は
G19を除いたこの環でも4である。上のfiber productでは最大idealの正冪が二因子の直和に
なるので、Hilbert–Samuel multiplicityは加法的である。従って

```
4=e(Ccore)=1+e(C_A),          e(C_A)=3.
```

平面hypersurfaceではe(C_A)=ord(g_A)なのでg_Aの最小全次数は3。
この証明は、旧abstract normal formのwを新しいphysical wと無条件に同一視したり、
g_Aを初めからsr²と置いたりしない。

## 4. 明示的なNoether parameter ell=5(s−r)

受理済みのactual conductor初期binary cubicと、その旧W,U,V tangent frameを使う。
旧29×22 frameから同じphysical G19 cutsのconstant Schurを引き、actual zero R0 sliceへ
投影した。旧UVはw速度零、旧Wはw速度非零であり、現在のz7の微分は旧UVで零となる。
全29行と七zero条件を確認した上で、旧(U,V)から現在の(r,s)への正確な行列は

```
[r]   [4/25  -4/25] [U]
[s] = [9/25   1/25] [V],       determinant=8/125.
```

従って

```
ell=−5r+5s=U+V  (一次部分).
```

旧V-axisではell=1、旧self tangent (U,V)=(1,3)ではell=4。
この二つがactual C_Aの接線であり、重みは単純・二重である。
受理済みの旧binary cubic `(12/125)U*(U−V/3)²` をell=0の方向(U,V)=(1,−1)へ
制限した係数は **64/375!=0**。actual g_Aの初期形との間の非零unit scaleは不要で、
ellがどちらのtangent componentにも消えないことが分かる。

ord(g_A)=3と合わせると、ellはC_AのNoether parameterでnonzerodivisor、

```
length C_A/(ell)=3.                            (3)
```

同値にC_AはC[[ell]]上finite free rank3。これはformal curveのWeierstrass theorem
またはDVR上の有限torsionfree性で従う。cubic Weierstrass係数を今回展開したわけではない。

## 5. T3の有限台torsionは一つのcolonで判定できる

(1)によりT3はC_A-module、また `(I0_CA,ell)` のradicalはR0の最大idealである。
従って

```
H0_mR0(T3)=H0_(ell)(T3).
```

ell-power torsionが非零なら、その適切なell冪を掛けることでell自体で消える非零元が得られる。
よってfinite-support torsionの **消失判定** は

```
(0:_(T3) ell)=0
 iff ((Q[3]:ell) intersect Jinf)/Q[3]=0          (4)
```

という一つのcolonだけでよい。colonは元のactual R0² submodulesのまま行う。
ell-kernelと元のm-socleの次元が同じという主張ではなく、その非零性・消失性の同値である。

このT3がfinite-support torsionを持てばfull横断Dもdepth0となり、cone depth19になる。
逆に(4)が零でも、T0,T1,T2の検査を省いてcone depth20とは言えない。
**今回(4)そのものは計算していない。**

## 6. 安全なbase changeと長さ上界36

I0_CAがT3を消すのでT3 tensor_R0 C_A=T3であるが、元の二行matrixを単にI0で
減らして像を取る操作は正しくない。Jinf tensor C_A→C_A²にはTorに由来するkernelがあり得る。
例えばR=C[[x]], J=(x), Q=(x²)ではJ/Q=Cは非零だが、J,Qのcolumnsをxで減らすと両方零になる。
必要なのは、T3自体のpresentationを得てからbase changeすること、または(4)を元のmodulesで
検査することである。96恒等式はこのTor問題を勝手に消していない。

T3はC_A上も十二generatorで生成されるため、(3)から

```
length(T3/ell*T3)<=12*3=36.                    (5)
```

これは閉fiberの長さ上界だけである。ell-torsionの指数や長さ、必要Taylor jet次数の上界ではない。
C[[t]]/(t^N)は一生成でmod t長さ1でもNが任意に大きいことが、その区別を示す。
したがって(5)を根拠にblindなjet次数増加を行わない。

## 保存・計算量

専用source `code/ambient_B18_T3_curve_annihilator_20261104_v63.py`。
`reports/AMBIENT_B18_T3_CURVE_ANNIHILATOR_20261104_v63/` に、96 pointer恒等式、
actual CA10×9 Jacobian、tangent transport、明示ell、actual g_Aの57unit input DAGを保存する。
一回0.997849秒、最大RSS141,576KiB。60秒/2GB内で、GB、jet、新点、core/source再生成は無い。
正確なT3 membership、colon、socle、全Dのnil/depthは次の独立した実計算である。
