# B18 actual係数DAG：J0二本と両六cubic値を構成、relative graph liftが残る

2026-11-04。要求した12個のsame-CA-trace-basis pairを**完成したとは報告しない**。
一回のbounded preflightで、actual CAのnu,J0二本、さらにCA/CBそれぞれの
全physical出力のw次数0..3係数を、有限DAGとunit-Jacobian代数方程式で構成した。
最後に残った障害は、CBの六cubic graphをCA targetへ評価するための、六つの
**具体的な1×10 odd-module membership**である。その実係数行列と右辺まで保存する。

## 1. 実際のzero15+base8 unit系を二次係数化する

共通physical zero baseは

```
R0=QQ[[x27,x108,x109,x112,x113,x117,x118,x119,x61]],
w=u140.
```

既存のcorrect core gaugeのsource未知数はzero group15とk0,k1,k3,...,k8の計23個。
これをzeta1,...,zeta23と書き、残るsource pを二次basisの元として

```
p²=s*p+t,        zeta_j=a_j+b_j*p
```

と置く。実際のzero gauge15式とphysical base9式、計24式の全てへ代入し、
basis(1,p)へ正確にreductionする。その48係数を零にする **実際の有限Hensel系**を保存した。
各式は保存済み150係数sourceと同じGL sectionを直接用いる。unknown係数の名前だけを
並べたrecipeではなく、四則演算DAGが48本の式そのものを指定する。

原点では全a,b,s,tが零。row順を偶24、奇24、column順を(a23,t,b23,s)にすると、
Jacobianはblock-lower-triangularで対角blockは二つとも

```
M24=[d_(gamma15,K8)(zero24 physical equations), [p²]source].
```

source p一次係数が零なので、偶blockにb,sの一次項は来ない。
奇blockのa,tへのcross項を零と仮定する必要もない。
新たにM24の全576 rational係数を保存してdeterminantを確認した。

|族|det M24|48次block determinant|
|---|---|---|
|CA|−9216/25|84934656/625|
|CB|576/25|331776/625|

これはcorrected23unit determinantとfold quadratic coefficientの積そのものである。
したがって形式的implicit function theoremにより、s,t,a_j,b_jはR0上の一意な
**algebraic power series**として定まる。既受理free quadratic normalizerのactual
basis(1,p)係数も同じ48式と原点条件を満たすため、一意性で一致する。
ここに未知のbranchやpの逆元を混ぜない。

```
tau=p−s/2,         nu=t+s²/4,
a+b*p=(a+b*s/2)+b*tau.
```

この式で同じ族内のtrace basisを全出力へ固定する。48式中のunit inversesは
全てorigin値を保存し、必要ならinverse-node uへu*v−1=0を加えてpolynomial systemへ戻せる。

## 2. positive五変数はall-ordersで線形に解ける

zero group G0はyと他の四座標を混ぜないblock diagonal matrixである。
positive group4方向は全てE_(i,y)、i!=y。そのため

```
G=G0+w*sum l_i E_(i,y),
G*G0^-1=Id+w*M,        M²=0,
```

となる。sourceのpositive変数はk2一本である。zero sourceはey wedge V(X)、
k2係数はyを含まない四座標上のbivector V2(X)であることを、**保存された全source項**から
検査した。sourceのweighted structureを後から仮定したのではない。

positive gauge4切断とphysical row140のw係数を使うと、l_iとk2/wの五未知数に対する
正確な5×5 linear systemがnormalizer内に得られる。その全pair matrixとGaussian traceを保存した。
両族でorigin pivot値は、row交換を含めて `(1,-1,-1,-1,1)`。
全pivotは二次algebra内のunitであり、p,s,tや未知のzero parameterを反転しない。

解いた後の全sourceはrank-one shearとして

```
(ey+w*l) wedge V(X−w*l*y) + w*k2hat*V2(X−w*l*y)
```

と書ける。V,V2はquadraticなので全w次数は正確に0..3だけである。
これは高次を打ち切ったsource jetではない。全150physical出力の四つのw係数を
各族の(1,tau) pairとして保存した。

## 3. Actual J0二本と六つのindividual cubic pair

共通positive odd physical rowsは **46,95**。その出力のw係数のtau係数が、
relative theoremに現れるactual J0二本になる。CAではscalar-DAG node IDs

```
j_A1=node17493, j_A2=node19106,
nu_A=node20086.
```

両jのorigin値は零である。これを係数函数の恒等零と混同しない。
CBの対応はnode38273,39886、nu_B=node40865。
zero odd I0八本も同時にactual scalar-DAG outputsとして構成した。

六つのphysical cubic rowsは `9,39,54,84,99,144`。
各族で

```
q_alpha/w³=A_i,alpha+B_i,alpha*tau_i
```

の全12 coefficient-functionsを構成した。
CA/CBの基底tau_iはそれぞれのactual quadratic normalizerから定義したものであり、
両者を同一視していない。

## 4. 要求されたrelative六pairに必要な不足データ

**individual physical cubic valueは、CB graphをCAのminimal targetへ代入した値ではない。**
両者のquadratic coverは同じphysical R0の上にあるが、tau_B=tau_Aという識別は無い。
二つのindividual pairを成分ごとに引くだけではactual交叉のgraph差にならない。

CBのzero oddを `q_i=A_Bi+f_Bi*tau_B`、positive oddを
`q_e=w*(A_Be+j_Be*tau_B)` と書く。H_B=R0+(I_B0+J_B0)tau_Bなので、各cubic odd係数について

```
B_B,alpha=sum_(i=1..8) h_alpha_i*f_Bi
               +sum_(e=1,2) h_alpha_e*j_Be       (★)
```

となるR0係数hが存在する。しかし、その **具体的なmembership coefficientsは未計算**。
既存のrelative packetは全87 graphの存在・weight boundを証明したが、そのscalar functionsや
この六つのmembership liftを保存していない。Noether fiberの旧source4は一変数上の有限jetで、
ここで必要な九変数の全次数membershipを代用しない。別々のabstract99モデルも代用しない。

今回のreader adapterは(★)の共通1×10行列の全node IDsと、六つの実際の右辺node IDsを
`relative-membership-target.json`へ保存する。これは未指定の未知係数を持つ入力ではなく、
**係数行列と右辺自体が指定された六つの有限algebraic-function module問題**である。
未知なのはそのsolution hである。

## 5. hが得られた場合のsame-CA-trace式を固定する

(★)を一つ解けば、CBのraw cubicを表す一つのactual graph liftは

```
Phi_Balpha=w³*A_Balpha
 +w³*sum_i h_alpha_i*(q_i−A_Bi)
 +w²*sum_e h_alpha_e*(q_e−w*A_Be).
```

これをCAの実際の同じphysical qへ代入してCA raw cubicを引けば、
CB−CAのnormalized coefficientは、**tau_A基底で**

```
a_alpha=A_Balpha−A_Aalpha
        +sum_i h_alpha_i*(A_Ai−A_Bi)
        +sum_e h_alpha_e*(A_Ae−A_Be),
b_alpha=−B_Aalpha+sum_i h_alpha_i*f_Ai
                    +sum_e h_alpha_e*j_Ae.
```

全てh以外の係数は今回のactual DAGで既に指定された。
要求された12組は

```
(nu_A*j_Ae*b_alpha, j_Ae*a_alpha), e=1,2; alpha=六rows
```

となる。しかしhを解いていないため、この12pairを完成出力とは呼ばない。
membership liftの選択差はCB kernel内の差であり、全CB kernelとlower graph relationsを
保持すればactual Q idealは変わらない。lower relationsを削除してこのtemplateだけを使わない。

## 6. 実行量と正確なfrontier

実際のsource preflightは一回、2.957198秒、最大RSS60,764KiB。
41,315 scalar-DAG nodes、certificate715,967bytesで、60秒/2GB/32MiB上限内。
ungauged R9のrank7計算、古いcoreの再生成、大きいGB、formal series展開、新しい点は無い。
正しいzero15+base8 selectorを使い、pの逆元も無い。

完成したもの：両48unit系、actual nu/I0/J0、positive5解、両族の全physical w0..3係数。
未完：六つの(★)のsolutions、same-CA基底の六relative graph pairs、12 topT3 pairsの所属、
四module socle test、全nilradical、depth19/20。

保存先 `reports/AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63/`。
主source `code/ambient_B18_actual_top_coefficient_preflight_20261104_v63.py`。
reader adapter `code/ambient_B18_relative_graph_membership_targets_20261104_v63.py` は保存済みDAGの
node参照を再整理するだけで、source calculationを再実行しない。
