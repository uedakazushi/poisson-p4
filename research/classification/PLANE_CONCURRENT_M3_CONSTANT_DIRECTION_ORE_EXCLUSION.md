# K0 constant primitive方向：全period/central/Wのactual因子排除

2026-11-08、PRIMARY G、producer純証明・独立review前。
新CAS、source probe、matrix、GB又はjetは無い。
全C_tan/C_transを閉じ、先行B_axis/B_otherと合わせ全constant-direction sourceを排除する。

## 1. 全linear-Casimir chartとPoisson vector field

Qはstrict binary sourceでnormal directionがconstant。
linear normalCasimir LでchartL1を取るとQ_L-column零なのでZ0。
全P=Q+∂c∧Wは

`P_chart=pi+∂c∧Y`, `L_Ypi=0`

で、Yは元全quadraticWから来るpolynomial chart vector fieldである。
actualDのchart値はpi wedge Y。
t0の全magnetic C_tan/C_transと、B_axis/B_otherは先行unimodular linear-Casimir lemmaで全空。
以下は残t非零を全て覆う。

## 2. C_tan t非零：polynomial transportがDをu³倍へ落とす

全t≠0でC2 shift map rank三なので、actualにC2零へ戻せる。
chartu1、r=v、x=b−arとする。このchangeはchart内のpolynomial automorphismで逆b=x+ar。
元πのliteral式から

`pi=∂a wedge V`, `V=∂r−t x∂x`。

全Poisson vector fieldYのtransverse componentを

`H=Y_x+t x Y_r`

とする。L_Ypiの(x,r)成分と残二成分から

`H_a=0`, `V(H)=−t H`。                               (1)

全YがpolynomialなのでHもpolynomial。
H=Σ h_j(r)x^jと書くとh'_j=t(j−1)h_j。
t非零ではpolynomial h_jが存在するのはj1、h1constantだけ。
従ってH=λx。
coordinate changeのJacobian一だからactualD chartはsignを除きλx。
quintic homogenizationは

`D=λ u³*(b u−a v)`。

λ零ならD零、非零ならproper factor。
全t非零/全Wを同時に排除する。
chart change自体をambient PGLだと仮定せず、得たactualDを元homogeneous coordinatesへ戻した。

## 3. C_trans：全parameterで一つの三角Darboux equation

chartv1、r=u、x=a/r−b、k=t+2とする。
一時r非零上で元πは

`pi=∂b wedge V`,
`V=r∂r−(k x+C2(r)/r)∂x`。

H=rY_x+(k x+C2/r)Y_rと置けば、全Poisson vector fieldから

`H_b=0`, `V(H)=(1−k)H`。                             (2)

元chart volume da db drとnew volume dx db drの比はrなので、
actualDの元chart値D0はsignを除きrH。
従って

`V(D0)=(2−k)D0`。                                     (3)

ここでpolynomial coordinate z=a−br=r xへ戻す。
(3)のderivationは元polynomial ring C[b,z,r]上で

`V=r∂r+(1−k)z∂z−C2(r)∂z`, V(b)=0。                 (4)

normalorder三のquinticDはbasedegree(a,b)≤二。
polynomialchange z=a−brはbasedegreeを保つのでD0のz/b総degreeも≤二。
さらに(2)はfixedx,rでD0のb非依存を与える（fixedz,rでも同じ）。
以下はこの必要条件全体だけを用い、元Yをchosen pointへ固定しない。

## 4. C2零：全kと全degree例外

C2零ならVはdiagonal。
z^i r^mのweightは(1−k)i+m、D0の必要weightは2−k。
z-degreei≤二だから全possible monomialsは

```
i0: r^(2−k),
i1: r z,
i2: r^k z²,
```

で、exponentがnonnegative integerの場合だけ存在する。
これで全複素kを扱い、integerだけをgenericと仮定しない。

u=rがDの因子でないためには、r-exponent零の項が必要。
上表からk0又はk2だけである。

k0ではD0はz²,r z,r²のspanに入り、元a,b,rでのordinary degreeは高々四。
quintichomogenizationはvで割れる。
k2ではr²z²が最高degree六のb²r⁴を持ち、他の二possible termsでは打ち消せない。
actualquinticのdegree≤五を戻すとその係数零。
残りはordinary degree≤三だから再びv因子。
他の全kではu因子。
zeroDもrank四既約対象外である。

## 5. 全nonzero central wallsはcommuting nilpotentで閉じる

先行exact shift mapによるC_transのnonzero residualは三wallsだけ。

| t | k | C2(r) | Vのnilpotent部分N |
|---|---:|---|---|
| −3 | −1 | r² | −r²∂z |
| −2 | 0 | r | −r∂z |
| −1 | 1 | 1 | −∂z |

nonzero coefficientはactual scaleで一へ戻した。zero residualは§4。
S=r∂r+(1−k)z∂zとすると全三場合で[S,N]=0。
Nは各有限polynomialを有限回で殺す。

polynomial D0をSのweight spacesへ分ける。
weightμ≠2−kでは(μ−(2−k))+Nがfinite nilpotent geometric seriesで可逆。
したがってD0はweight2−kに入り、ND0=0。
domain性から∂zD0=0。
(2)のb非依存も戻すとD0はconstant倍のr^(2−k)だけ。
指数は三場合で3,2,1だから、全てu因子。
resonanceをgenericfieldの分母で除外したのではなく、全wallのwholeDを直接扱った。

## 6. 結論

C_tan/C_transの全period、全centralzero/nonzero、全quadraticWは、
actualD零又はnormal hyperplane proper factorを持つ。
先行unimodular lemmaのB_axis/B_otherと合わせ、
**K0の全constant primitive direction源は今回の既約quinticを持たない**。
残るK0はprimitive invertible linear directionのSS/Jordanだけ。
ここで得たpolynomial chart normalformsをambient actual軌道分類とは呼ばない。
SS/Jordan、wholeconcurrentG、E又はcomponent数の結論はまだ残す。
