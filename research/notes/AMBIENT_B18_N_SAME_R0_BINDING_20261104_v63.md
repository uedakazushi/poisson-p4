# Actual N tailを同じphysical R0上のR8+m8 tauへ結合する

2026-11-04。個別Ac7型から座標を推測せず、既存のN normalization sourceとliteral W1 cut mapを
結合して、実際のR0九関数からCohen baseと残るodd座標を選んだ。
Jinfのsyzygy、D4のextension relation、T3の零性はまだ未計算。

## 1. 同じN sourceと実際の20切断

凍結済みN normalization sourceは(p,s,k0,...,k8)の11変数で、中心は
p=s=0,k0=1,k5=1/3、他k=0。全61出力は既知の厳密な多項式である。
これを旧49表示へ、ss=s,aa=p²,bb=p³,cc=spと実q45のまま代入した。
前packetの同じrational group17を加え、physical20切断とR0九関数を読み戻した。
groupのbranchや別のbase座標を暗に選んでいない。

normalization28変数に対する20切断のunit minorは **-96**。
従ってそのnormalization tailはregular8次元である。
R0九関数の一次微分rankは **7**。選んだ七つの実際のR0関数はphysical rows

```
108,109,112,113,117,118,119.
```

これらをX0,...,X6と書く。20切断、この七関数、元sourceのpを合わせた28系の
Jacobian determinantは **23040**。
従って(X0,...,X6,kappa=p)はactual regular normalizationの形式座標である。

旧normalization quotientの自由3変数(k0,k2,k5)とgroup17の計20変数に切断を制限した
minorも **-480** である。これはnonnormal smooth factorが切断で失われることの実確認である。
以下のimage idealの証明は、さらに直接のodd unit matrixも用いる。

## 2. Foldと九番目の関数を二次の必要成分だけで測る

20切断と七Xを固定した28-unit系で、kappa²とkappa*Xiを保つ二次quotientを計算した。
保つmonomialは1、八つの一次項、七つのXi*kappa、kappa²だけ。
Xi*Xjと三次以上を落とす明示的なquotient algebraなので、必要成分への正確な計算である。
任意の高次jetや点探索はしていない。

残る二つのphysical R0関数はrows27と61であり、kappa²係数はそれぞれ

```
rho27: -3/4,        rho61: -3/2.
```

従って八番目のCohen parameterには **V=rho27** を取れる。
actual R8=QQ[[X0,...,X6,V]]は元R0の八つの座標そのものであり、normalizationは
その上finite free rank2となる。これはregular sourceのkappa方向で、Vがunit*kappa²から
始まることに対するWeierstrass theoremである。

九番目rho61の、Vの二次寄与を差し引いたoddの一次係数は

```
b_linear = -(1/2)X2 + 2X4
         = -(1/2)rho112 + 2rho117.
```

実際のtracefree basisへの変更を含めてもこのmixed二次の式は変わらない。
この時点でbがregular parameter候補であることが分かるが、以下で全有限inputも構成した。

## 3. 純粋な56-unit係数系がa,b,nuを指定する

source pをquadratic basisとし、他の27個のsource変数をAi+Bi*p、
p²=S*p+Tと置く。Biの中心は上のunit28系から得たkappa-kernel係数を用いる。
20切断と七つのX、およびrho27=Vの計28式をbasis(1,p)へ正確に戻すと、
56係数未知数(A27,B27,S,T)への56式となる。
原点Jacobian determinantは **298598400**。
全系・全Jacobianをliteral scalar DAGとして保存した。

全分母は原点で単位であり、unit Henselの一意解がR8上のalgebraic power seriesを定める。
sourceのrank2 Weierstrass表示と同じ中心を持つので、この系はその係数を一意に回復する。

```
tau = p-S/2,        nu=T+S²/4,
rho61 = a(R8)+b(R8)*tau.
```

nu,a,bは証明書に実際のnodeとして指定した。
fold係数-3/4から、dnu/dV=-4/3が非零であり、normalizationはregularである。
56系のimplicit ADでもbの全八parameter微分を計算し、上のmixed二次式と一致した。
V方向のb一次係数は0である。

## 4. Actual targetのodd idealは正確にm8

全150の実physical target関数を同じquadratic系で評価した。
tracefree basisで全odd係数の原点は0であり、actual targetはR8+m8*tauに入る。
以下八つのphysical rowsのodd係数を順にz1,...,z8とする。

```
61,18,19,28,29,45,46,60.
```

これらのR8八parameterへのimplicit Jacobian determinantは **-4/9**。
従ってactual odd係数ideal IはI⊂m8、I+m8²=m8を満たす。
NakayamaからI=m8が全次数で従う。
targetはCohen base R8自身を含むので、各targetのeven成分をR8関数として差し引ける。
またnormalizationはfinite R8-moduleだから、image algebraも有限moduleで完備閉である。
従ってactual completed N imageは正確に

```
A_N = R8 + m8*tau ⊂ R8[tau]/(tau²-nu).
```

これは個別Ac7という型からの推測ではなく、全actual targetと同じR0上のunit8行列による証明。
normalization sourceのproperな唯一の幾何fiberと同source chartのunit条件を使うため、
別枝を落としたimageではない。

## 5. 同じR0の九番目はv=z1*tau

選んだR8は元physical R0の八座標なので、元R0はその上rho61を加えたregular9変数環。
rho61のeven部分a(R8)を引く変更は、R0自体のunit formal automorphismである。

```
v = rho61-a(R8),       R0=R8[[v]],
v maps to b*tau = z1*tau in A_N.
```

z1=bに続く上の八odd関数はR8のunit座標系である。
全56係数系とz_i=b_i(X,V)という追加8 unit equationsを合わせて、
nuとaもこの新しいz座標上の有限implicit関数として指定できる。
元R0を別のabstract baseへ無断に置き換えたわけではない。

したがってNの同R0 moduleには、generator1,z_i*tau (i=1..8)、Koszul28と

```
v*1 = z1*tau,
v*(zi*tau) = nu*z1*zi
```

という具体的なsmall presentationが適用できる。
最初のunit relationでz1*tauを消すと8 generator、36関係となる。
このpacketはその全係数とactual target generatorのmapを指定するが、
36行列・次syzygy・D4へのextension liftsのliteral構築は次の独立した手順とする。

旧nil curveのphysical parameterはrho117だった。
bのこの方向の微分は2なので、新しいz1も同nil curveのunit parameterとなる。
これはD4 extensionを一変数上で調べる際の同R0 compatibilityを与えるが、
extension coefficientやF_oldのH_A二行liftを求めたことではない。

## 6. 範囲と保存

tailのw=1操作はfinite-w polynomial algebraを先に局所化してから行う。
`AMBIENT_B18_W1_POLYNOMIAL_LAURENT_SCOPE_20261104_v63`の訂正を適用し、
w=0完備環そのものへw=1を代入しない。

Source `code/ambient_B18_N_same_R0_binding_20261104_v63.py`。
`reports/AMBIENT_B18_N_SAME_R0_BINDING_20261104_v63/actual-N-R0-binding-DAG.json`に
全55357 scalar node、Nsource28/実cut20、R0九関数、unit28frame、必要な二次係数、
全56-unit系、全150odd微分、unit8座標変更を保存した。
一回0.938951秒、最大RSS83404KiB、証明書982280bytes。
GB、新点、三次以上のsource jet、opaque branch pullback、p逆元は無い。
T3の零性、全B18原点のdepth/nil、D4 extensionの分裂は主張しない。
