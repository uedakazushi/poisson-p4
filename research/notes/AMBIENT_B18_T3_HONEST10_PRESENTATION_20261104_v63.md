# Actual D3/D4の正確な10行presentation：Torを捨てない入力

2026-11-04。全186入力の二行matrixをactual CAの正確な10-generator moduleへ持ち上げた。
T3のrelation module、零性、depth19/20は依然未計算である。

## 1. Hの全kernelは36本のKoszul関係

R0はactual共通physical九変数regular完備局所環であり、actual CAの九座標を

```
Y=(f18,f19,f28,f29,f60,f67,f72,j46,j95)=(z7,r,s)
```

とする。determinant128/3が受理済みなのでm_R0=(Y1,...,Y9)。
次数d>=1のCA component moduleは

```
H=R0 + m_R0*tau_A.
```

10個のfree generator E0,...,E9を
`E0 -> (1,0), Ei -> (0,Yi)`と元のtrace pairへ送る。
normalizerがbasis(1,tau_A)でfree2なので、kernelはeven成分0、oddのY-syzygyだけである。
Yはregular sequenceだから、その全first syzygyは36本

```
-Yj*Ei + Yi*Ej,    1<=i<j<=9
```

で生成される。これらを保存した。したがって元の二行matrixとは異なり、
R0^10からHへのkernelも全て明示された正確なmodule presentationである。

## 2. 古いg商をrestrictionして0/0にしない

全186 scalar DAGにはcoordinate商とg商がある。これらをprefix inverseへ機械的に代入し、
分子・分母とも恒等零となったものを再除算する方法は採らない。

代わりに、CA自身の元のpure-unit source functionだけを座標分解した。
必要なものはf_CA63と87 physical graphのCA odd coefficientの計88個。
全て原点0であり、既存CA48系とF9からなるgenuine unit57系で定義される。
CA F9のprefixを一つずつ残す8個の57系を作る。
九つ全てを残したinverseは元のsourceとpointed uniquenessで一致するので、新しい系は不要。

各関数bについて

```
h_i=(b(Y1,...,Yi,0,...,0)-b(Y1,...,Y_(i-1),0,...,0))/Yi
```

を792個のregular coordinate quotientとして保存した。
二つのprefix系はYi=0で一致するから分子はYiで割れ、algebraic power seriesのままである。
全式がtelescopingし、`b=sum Yi*h_i`を全次数で与える。
各57系のdeterminantは−3623878656/625。
新しい商の原点値はactual source ADとCA F9逆Jacobianから計算し、全保存した。
数値jetの零からmembershipを推論したものではない。

## 3. 元のfactorized DAGに沿う正確な持ち上げ

CA zero7とj46,j95には座標basisのmembership vectorがあり、上の88関数には
telescoping vectorがある。これらから次の厳密な規則でodd membershipを伝播した。

- 和・負号：membership vectorの和・負号。
- 積：一方のfactorに既知membershipがあれば、他方のscalar functionを掛ける。
- 元のDAGが`-(-b)`をbへ簡約していた場合：既知関数の反対符号も登録する。

CB側coefficientの古いr/g関数は、この積のscalar cofactorとして元のphysical base上に
保つ。その関数をrestriction・再除算・未知微分0へ置換していない。
全CB kernel odd termはCA odd factorを含み、各graph odd termは
`sum h_CB,i * CAodd_i - CAgraphOdd`であるため、この規則が全てに成功した。

保存済みQ[3]2034列、top12列、Jinf2046列の全てのodd coefficient bに、
9成分vector vで`b=sum Yi*v_i`という全次数の恒等式を与えた。
従って元のtrace pair(a,b)のliftは正確に`(a,v1,...,v9)`である。
全propagation ruleを関数nodeごとに記録した。

## 4. 得られた実際の有限入力

```
D3 = coker( R0^2070 -> R0^10 ),
 2070列 = 36 Koszul + 2034 actual Q[3] lifts.

D4 = coker( R0^2082 -> R0^10 ),
 2082列 = 上の2070列 + 12 actual top lifts.
```

degree4ではJinfに安定し、Jinf=Q[3]+top12であることは元のdegree別column分解と一致する。
T3はD3内の12 lifted top columnsの像であり、正確に`ker(D3 -> D4)`となる。
これは新しいactual constructive inputであり、単に抽象的な表示の存在を述べたものではない。

ただしT3そのもののpresentationは、まだ求めていない。
そのためには上の2070関係に対するtop12のrelation module、またはJinfの正確なsyzygyを
計算する必要がある。D3のmatrixをm²で減らすだけではT3のpresentationにはならない。
既受理のm²*T3=0も、この最後のkernel計算を自動的に消すものではない。
全nilradical、T3非零、depth、colon、positive membershipは今回主張しない。

## 保存・実測

Source `code/ambient_B18_T3_honest10_presentation_20261104_v63.py`。
`reports/AMBIENT_B18_T3_HONEST10_PRESENTATION_20261104_v63/honest10-presentation-DAG.json`に
全486437 scalar node、8×57unit equations、792regular quotients、全membership規則、
36Koszul列と全lifted matricesを保存した。証明書10065275bytes。
一回1.723826秒、最大RSS255148KiB。
最初の構築試行では既存DAGの二重負号簡約で既知関数のchildを認識できず停止し、
厳密な反対符号登録を加えてから成功した。sourceや計算次数を変更した再探索ではない。
GB、jet、新点、source/core再生成、p逆元、old g商のfunctorial pullbackは無い。
