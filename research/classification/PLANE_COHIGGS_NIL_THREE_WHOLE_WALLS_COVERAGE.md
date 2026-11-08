# Co-Higgs rank M三nil：残三whole壁の点ごとの被覆

2026-11-08、PRIMARY E、log_manuscript。producer全源証明、独立受理は別義務。
1054主稿／frontier／global／Git／旧bytesは不変。
新component数の増分零。全rank零plane又は全Eの完成は主張しない。

## 1. 全sourceと唯一の有限batch

[有限入口](PLANE_COHIGGS_NIL_THREE_WHOLE_WALLS_ENTRANCE.md)の元全三壁
tau=−3/2,1/2,3/2を扱う。
全rank三Lax正規化と実逆、全P1三十、FULL九P2、全元350式を保持する。
保存済みSの同じminorだけでK=S A+N Bを全30行で分解し、
旧D2の同じminor inverseでP2を復元した。
新M1、rank、minor又はpointは作っていない。

tau半ではS rank五、全kerD1七なので二つのNを保持し、
さらにD2核一を自由係数rhoとして完全に保持する。
他二壁のS rank六、D2単射もそのまま使う。
全六actual shearのpolynomial正逆、保存SとのEuler商等式、
全parameterの150係数push/inverse roundtripを確認した。

[inline job](jobs/plane_cohiggs_nil_three_whole_walls.json) のgate PASS後、
180秒／3GiB契約でwhole三源を一回だけ実行した。
実測0.208868秒、最大RSS14,760KiB。
[不変report](../../reports/CLASSIFICATION_PLANE_COHIGGS_NIL_THREE_WHOLE_WALLS_20261108.json)は
実行source全文、全N、全九P2、全350 residualと全126 actual Pfaffianを保存する。
producer／CASをimportせず、既独立Fraction sparse primitivesだけを再利用した。

最初の有限reportのautomatic criterionは、tau半の一faceに
UNRESOLVED_exact_whole_parameter_faceを正直に残し、
all_three_whole_sources_closed=falseを保持している。
これはcoordinate factorだけの有限診断であり、§4の非coordinate cubic因子はその後の純展開証明である。
旧falseをtrueへ書き換えず、その別証明を接続して全sourceを閉じる。

## 2. tau負三半：全一parameter源はc二乗で割れる

残normalをlambdaとし、唯一のNの非零括弧は

~~~
Nab=aw, Nbc=−cw/2, Nbz=zw/2.
~~~

Nは共通∂b方向とb非依存係数を持ち、strict Jac(N)=0。
旧D2単射から全P2零で、全模型はP0+lambda N。
元actual五次は

~~~
D=−c²(3q+lambda w³)/2,   q=az²−2bzw−cw².
~~~

qのa z²係数は一なので第二因子も全lambdaで非零非定数。
従って全lambdaの指定既約quintic源は空となる。
全sourceのshear逆／Lax normalizer逆を戻しても同じである。

## 3. tau三半：全二parameterに第二m三Poisson plane

全残Nは次の二方向で、係数をlambda,muとする。

~~~
N1ab=cw, N1aw=2zw, N1bw=z²;
N2ab=−9bw, N2ac=bz−9cw/2.
~~~

FULL九P2の必要唯一解は

~~~
P2ab=−2lambda² z²+9lambda mu zw,
P2ac=−lambda mu z²,
~~~

他の全P2零である。全350は全lambda,mu上で恒等零で、
全parameterがactual Poisson tensorを与える。

S'=(c=z=0)は、全c-row／z-rowを代入すると全parameterでPoisson plane。
そのrestriction cubicは

~~~
R'=4b²w−9mu bw².
~~~

固定b²w係数四が非零なので全parameterで一般rank二。
actual五次のnormal位数は正確に三で、最低symbolは

~~~
G3=−3c³w²/2−3bc²zw−9mu bcz²w/2−mu b²z³.
~~~

c³w²係数負三半はparameterに依存せず非零。
lambdaはこの最低symbolを変えない。
mu非零のreduced concurrent cubicとmu零の2L+Mの両境界を含み、
既受理の全rank二plane／m三定理を全sourceへ適用できる。
よって全被約既約quintic点は既知閉成分の有限和に入る。
zero境界をgeneric Rへ押し込んだのではない。

## 4. tau半：元uv行で二whole faceに尽くす

全二残Nの係数をlambda,mu、旧D2核の自由係数をrhoとする。
非零括弧は

~~~
N1ac=2bw, N1az=w², N1bc=bz−cw/2;
N2ab=bz, N2aw=−3zw, N2bw=−z²/2.
~~~

FULL九P2は正確に

~~~
P2ab=mu²z²+(2lambda mu−rho)w²,
P2ac=(2rho+2lambda mu/3)zw,
P2bc=rho z².
~~~

rhoを零又はgenericに固定していない。
元projective Jacobiのrow51は正確に−lambda muで、
残八非零rowもlambda muの倍数である。
従って元全Poisson supportはlambda零又はmu零の二閉faceに尽くされる。
両faceで全350が零であることも同じ全係数から分かる。

### lambda零：第二planeの全m三帰属

S'=(c=z=0)は全mu,rhoでPoisson plane。
restriction cubicと最低normal symbolは

~~~
R'=2b²w−rho w³,
G3=−w(c²−2rho z²)(bz+cw/2).
~~~

それぞれb²w係数二、c³w²係数負半が固定非零なので、
一般restriction rank二かつnormal位数三を全parameterで保つ。
rho非零／零で零cubicとnormal根型が変わる全境界も、
既全m三pointwise定理が含んでいる。
よって全mu,rhoの指定実既約quintic源は既知閉和に入る。

### mu零：全actual五次に非coordinate cubic因子

残る全lambda,rho模型について、元actual五次は手で展開して正確に

~~~
D=(q−lambda w³)*
  (c²/2+lambda c w/2−lambda b z−rho z²),
q=az²−2bzw−cw².
~~~

第一因子のa z²係数一、第二因子のc²係数半が全parameterで非零。
従って両方が常に正次数の非零多項式であり、
全lambda,rhoの指定被約既約quintic源は空である。
有限reportがcoordinate-factor診断で認識しなかったのはこのcubic因子だけで、
新anchor／matrix／jet／GBを足したのではない。

この恒等式は、旧全Pfにmu零を代入した全多項式と二因子の積を
標準Fraction疎多項式でliteralに比較する補助readbackにも保存する。
その一つの乗算恒等式は新operator又は別のwhole-source実行ではない。

## 5. 元全sourceへの戻しとrank三全体

全三壁のK=S A+N Bは完全核の全係数を表し、
finite shearは全P2を同時に動かす実global linear changeである。
したがって§2–4の全factor／第二plane結論を全元sourceへ戻せる。
known unionはPGL不変かつ閉であり、個々の成分番号の識別は不要。

既受理[全m三pointwise theorem](PLANE_RANKTWO_M3_WHOLE_POINTWISE_COVERAGE.md)を使うのは、
実tensorの第二Poisson planeとrank／multiplicityを全sourceで明示した後だけである。
五次式の一致からtensorの親を推測したものではない。

[semisimple全tauとnil tau一の純証明](PLANE_COHIGGS_FULLRANK_M2_SEMISIMPLE_AND_NIL_ONE.md)、
generic／負半の全shear-exact空性、
[tau零の全rank二定理](PLANE_COHIGGS_FULLRANK_HIGH_M_WHOLE_RANK_TWO.md)と合わせると、
rank零planeのO(1) co-Higgsでrank M三の全multiplicity源は、
指定実既約quinticの下で既知閉成分の有限和又は空へ全点被覆される。
rank M二M(p)非零のm二壁、rank M一以下、scalar、zero-symbol等は別の残余である。
このaggregateの独立受理／採録と主稿統合はrootへ委ねる。
