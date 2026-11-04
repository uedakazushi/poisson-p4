# Actual B18 full相対入力：全99+87関係と四つの有限module行列を構成

2026-11-04。前のactual source DAGと単元9座標を一般化し、全てのCB kernel関係と
全87 physical graph差を、**同じCAのtrace基底**で有限algebraic circuitsへした。
これはfull Dの係数入力の構成であり、全nilradical・primary分解・socle・深さ19/20を
解いたという主張ではない。これまでのblocked stageとtop12 stageは変更しない。

## 1. 実座標とgの非零性

物理base R0の九変数、actual CA/CBの48unit Hensel systems、source(1,p)からの
trace basis `tau_i=p_i−s_i/2` を前二packetのまま使う。
CBの選択座標は

```
(z1,...,z7,r,s)
 =(f18,f19,f28,f29,f60,f67,f72,j46,j95).
```

R0におけるJacobian determinantは1/96。残るzero odd coefficientはf63。
unit inverseで最初の七座標を零にして

```
g(r,s)=f63(X(z=0,r,s))
```

を定義する。七つのcoordinate telescopingで実際に

```
f63=g+sum_(i=1..7) zi*ei
```

という係数eiを構成した。gが零ならactual odd ideal
`I_CB=I_CB,0+w J_CB,0` は七つのziとwr,wsの九本で生成される。
これは独立受理済みの**actual mu(I_CB)=10**に反する。従ってgは非零で、
g(0,0)=0、gはw,zに依存しない。

従って、実際のgenerator changeによって

```
I_CB=(z1,...,z7,g(r,s),wr,ws).                  (1)
```

となる。元のf63 targetからsum ei*(zero odd targets)を引くunipotentな変換であり、
abstract normal formの座標を相対的に同一視してはいない。

## 2. 四つのzero graphのregular divisionに循環は無い

独立なcanonical ambient126座標を、既存のcommon minimal39 physical座標から延長した
凍結済みreaderの87行を使う。zero-weight rowsは `64,68,73,74`。
CBの各zero physical出力をa_alpha+b_alpha*tau_Bとする。
そのb_alphaがI_CB,0に属することは、未構成のrelative graphを仮定するのではなく、
**actual imageのminimal algebra generation**から分かる。

すなわち、共通minimal20のphysical functionsはcomponent完成環のcotangent basisであり、
形式的Nakayama/Cohen presentationにより全完成像algebraを生成する。
nonnegative gradingによりzero部分はbase R0と八つのzero odd physical functionsだけで生成される。
各々のodd coefficientはf_iであり、その任意の積のodd coefficientはideal(f_i)内。
formal極限でもNoetherian idealは閉じているため、全zero出力のb_alphaはI_CB,0内にある。

z=0へ制限すると

```
bbar_alpha(r,s)=b_alpha(X(0,r,s)) in (g).
```

よって

```
h8_alpha=RegDiv(bbar_alpha,g)
```

は正則algebraic power seriesとして一意に定まる。numeratorとdenominatorは共に今回の
有限DAGで指定され、未知係数を立てたmembership recipeではない。
g!=0と上の既証membershipが正則性を保証する。pや任意の未知factorの逆元ではない。

七coordinate方向のtelescopingをb_alphaに適用した係数をh_alpha,i^zとすると

```
b_alpha=h8_alpha*f63
        +sum zi*(h_alpha,i^z−h8_alpha*ei).
```

これで四つ全てのzero graphのactual odd membership liftが得られる。
この四RegDiv(g)の**origin値だけは未評価**だが、関数自体は上式で完全に指定されている。
それらを単元として使用しておらず、未評価の値を零と置くこともない。

## 3. 全83 positive graphも同じ九座標で持ち上がる

weight1/2/3のphysical graphは53/24/6本である。
各normalized出力 `q_alpha/w^d=A_Balpha+B_Balpha*tau_B` のBは原点で零で、
actual `I_CB,0+J_CB,0=m_R0` に属する。
前packetの九coordinate telescopingを全83本へ適用する。
zero7を消してr,sを残すinverse instancesと合わせ、必要なunit57 systemsは計16個。
coordinate regular divisions784個、既証membershipによるRegDiv(g)四個を保存した。

これらのliftを同じphysical minimal20 functionsに書き、CAの実際の出力へ代入する。
各族のtauを同一視せず、R0上のscalar coefficientsだけを共通physical baseで評価する。
全87個のnormalized CB−CA graph pairを **(1,tau_A)** 基底で保存した。
degree分布は正確に4/53/24/6であり、一つも削除していない。

## 4. CBの全44 first syzygiesを有限に構成する

g(0,0)=0なので、もう一つの有限coordinate splitにより

```
g=r*A+s*B,
A=(g(r,s)−g(0,s))/r, B=g(0,s)/s
```

を正則algebraic functionsとして得る。gのfactorizationやsocle計算は不要。
R0[[w]]のregular coordinates(z7,r,s,w)を使う。
core ideal (g,wr,ws) はheight2で、二つの全syzygiesは

```
(w,−A,−B),   (0,s,−r).
```

直接にも証明できる。`a*g+b*wr+c*ws=0` をwで割った商環で見るとg!=0なのでw|a。
a=w*a'としてg=rA+sBを入れれば、rとsのregular-sequence syzygyから上の二本で尽くされる。
これはHilbert–Burch matrix

```
[[ w, 0 ],
 [−A, s ],
 [−B,−r ]]
```

のminors(g,−wr,ws)でもある。

(1)の七つのziを戻すと、mod(z7)で上のcore syzygiesを引き、ziを含むcore係数を
cross relationsで移し、最後にregular z7のKoszulを使えば全first syzygyが尽くされる。
従って完全な生成系は

```
21 regular Koszul +21 zero/core cross +2 core HB =44.
```

old f63へ戻すunipotent generator changeをactual targetへ施す。
十個のodd targetsをZ_i、対応scalar coefficientsをb_iとすれば、さらに55本

```
Z_i Z_j−nu_B*b_i*b_j (i<=j)
```

を加えた **全99本がCB image kernel全体を生成**する。
証明はproductsで全式をR+sum R Z_iへ還元し、basis(1,tau_B)で残る条件が44個で
生成されるfirst syzygyであることによる。99本の最小性は主張しない。

## 5. 全186 Q入力をsame CA trace basisで評価する

各CB odd targetを、CA raw physical出力からCBのeven trace partを引いた実際のpairへ送る。
f63のunipotent changeもこの同じ入力へ施す。
product relationsの右辺には **nu_B**、pair multiplicationには **nu_A** を使う。
両者を取り違えず、tau_A²=nu_Aのもとで全99関係をR0² pairへ評価した。

87 graphとの差を合わせ、actual D=A_CA/QのQは、指定された全186 homogeneous generatorsで尽くされる。
それぞれはw^dを除いたsame-CA-trace pairとして保存し、degree分布は

```
CB kernel :64,32,3,0
graph     : 4,53,24,6
total     :68,85,27,6  (degrees0,1,2,3).
```

完成像ringそのものは既知のactual
`A_CA=R0[[w]]+(I_A0+wJ_A0)R0[[w]]*tau_A` を用いる。
二normalizerのtensor productをactual intersectionと見なしてはいない。

## 6. 四module criterionの実際の有限行列

actual CAのC0 generatorsは(1,0)と八つの(0,f_Ai)、H generatorsはこれに
二つの(0,j_Ae)を加えたもの。これらを全186pairへ掛け、R0²内の実際の生成行列を保存した。
以前の四module criterionの記号で

|module|保存した生成column数|
|---|---:|
|Q[0]|612|
|Q[1]|1513|
|Q[2]|1926|
|Q[3]|2034|
|Jinf|2046|

全て有限algebraic-function DAG entriesであり、未指定coefficientを持つsymbolic matrixではない。
Jinf/Q[3]の追加12columnsも指定した。生成系の最小性やmatrix rankを主張しない。

これで次の実際の問題は、指定済みの四module socle/colonと上記12columnsの所属である。
非coordinate regular quotientを普通のunit inverseとして扱ったり、pointed Hensel branchを
失った未局所化polynomial GBへ移したりしてはならない。
現時点で所属・colon・socleは計算しておらず、depth19/20・nilradicalは未決である。

## 保存

専用source `code/ambient_B18_full_relative_input_20261104_v63.py`。
`reports/AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63/full-relative186-DAG.json` に
全16unit57 systems、regular quotients、99 kernel、87 graphs、186pairs、module行列を保存。
一回1.185516秒、339,669nodes、6,769,682bytes、最大RSS169,772KiB。
60秒/2GB/32MiB内で、source/core再生成、新点、GB、series truncationはない。
四RegDiv(g)の正則性は§2の既証minimal algebra generationに依存することをreceiptにも明記した。
