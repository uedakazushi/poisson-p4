# 非被約R・m二moving doubleの独立全P0空性監査

2026-11-08、PRIMARY G、`/root/log_manuscript`。
rootの `PLANE_NONREDUCED_M2_MOVING_DOUBLE_EXCLUSION.md`
（SHA256 `7dec4418bb8dcd8df242602940c01230dabf35cf61abbfd2dde31705b6f9c58e`）を
別作者として読解し、全support、nonlog trace、vertex条件及び全actualP0係数を独立に確認した。**PASS**。
新CAS、matrix、probe又はjetは実行していない。

## 1. 全leadingとseedへの正当な必要条件

対象は実D被約既約quintic、restriction rank二plane、f=a³又はa²b、normal位数二、
非零double G2、primitive root degree一の全源である。
UFD/二gradingで `G2=gamma*h*(Au+Bv)²`、gamma非零、h一次、A,Bはindependent一次、gcd一。
hはR支持なのでa又はbに限り、q=f/hは二次。
weight零P0自体がPoissonで、actualD(P0)は全G2そのもの。
Q=Au+Bvはrank四、vertexはp=(A=B=0)であり、hはQの因子ではない。

既受理 `DOUBLE_QUADRIC_SUPPORT.md` §4の「rank四Q、単純Hがvertexを避ける全源は空」は、
このP0に仮定を満たして適用できる。ゆえにp∈h。
P0の必要条件を使っているだけで、seedの既成分帰属をunknown higherへ延長していない。

## 2. 唯一の局所単極比較と非被約成分の保持

qがhを含めばp∈qは自動。
残るf=a²b,h=b,q=a²で、もしa(p)非零ならp近傍のfはunit×bである。
**この近傍だけ**でbetaは高々単極で、離れたa²成分の極を単極と仮定しない。
tracefree Gammaとroot rowは

```
d ell+ell Gamma=(dlog q/2)ell.
```

p∈bよりbはA,Bのspan内なのでconstant normal右GL2とlocal base coordinatesでell=(b,t)と出来る。
Gamma=(U db+V dt)/b、regular部分もU,Vに吸収する。
bを掛けてb零上へ制限したroot式は、両residue matrixの第二rowを零にする。
trace零と合わせU,Vはstrict upper triangular。
flatnessのb^-2 dbdt係数はVのresidueと[U,V]の符号付き差で、後者零からV|b0零。
よってdt成分Gamma_dtはregularであり、pのdt root式は左辺(0,1)、右辺及びell Gamma_dtは零。
矛盾するからp∈qも全源で必須である。
irregular traceはtracefree式で消えており、その係数を零に設定したのではない。

## 3. 固定traceの全connection一意性

T=tr betaを固定すると、同じprimitive ellを持つ二connectionの差は

```
delta beta=omega*[[-AB,-B²],[A²,AB]].
```

matrix entriesのgcd一とf delta betaのpolynomial性から、f omegaもpolynomial one-form。
degree比較でそのcoefficientsは定数、Euler beta=IからiE(f omega)=0なので全て零。
この議論はf被約を仮定しない。

p∈qはq∈(A,B)と同値なので、linear auxiliary rowを補いdetM=qと出来る。
`theta=(T+dlogq)/2`、`beta=theta I−M^-1dM` は元root、trace T、Euler Iを持つ。
closed TとMaurer--Cartan式からrow conventionのflatness d beta−beta∧beta零。
f thetaはpolynomial、f M^-1dM=h adj(M)dMもpolynomialなので、全global sourceに属する。
ゆえに一意性から元betaそのものである。
auxiliary rowはprimitiveでなくてもよく、actual rootをmixするleftGL2を同値として使わない。
Tの全irregular exact部分を保持した全traceごとの一意性である。

## 4. 全q=abとq=a²のactualP0矛盾

q=abの場合はf=a²b,h=aであり、p∈aからactual normal右GL2でell=(a,t)、t independent aへ戻せる。
completion `M=[[a,t],[0,b]]` のbeta21は零。
canonical tracefree mixed gaugeのnormal-v columnはv倍となるから、
actualG2のu²係数は、全normalKを残しても正確に3f k20だけである。
target gamma*h*a²との比較はb|aを要求して矛盾する。
tがbと比例する場合、pがa∩b又はaの他点の場合も同じ非整除であり、位置をgenericにしていない。

q=a²ではf=a³,h=a、又はf=a²b,h=bの両方を含む。
p∈aから同じell=(a,t)として `M=[[a,t],[0,a]]` を取ると、
betaのdiagonalは同一、off-diagonalは上側だけである。
したがって全raw sourceは

```
P_base,u=uX+vZ, P_base,v=vX,
P_uv=k20u²+k11uv+k02v².
```

actual uv係数のmixed determinantはX∧Xで零、二normal Euler寄与はX(f)の差で零、Zはv²にのみ寄与する。
従って **actualG2[uv]=3f k11**。
targetの係数2gamma*h*a*tと比較すれば `3a k11=2gamma*t` となり、a,tのindependenceに反する。
factor三はP_fのE(f)=3fから、factor二は実square rootの交差項から生じる。
normalKの他二係数やirregular thetaでこの非整除を解消できない。

以上で全moving degree一double rootのP0が空、従って全higher源も空。
定数root degree零、SF、fの他型、normal位数三及び全G/Eはこのreceiptに含めない。
全P1三十・P2九を狭めておらず、P0の必要係数だけで矛盾した。
