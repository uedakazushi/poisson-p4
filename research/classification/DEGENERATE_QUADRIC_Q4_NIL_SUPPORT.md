# rank四の二次因子：nil(3,1)・nil(2,2)五源の全支持還元

2026-11-07。課題 `E_reduced_nonnormal`。
対象は `DEGENERATE_QUADRIC_CUBIC_COFACTOR_REDUCTION.md` のrank四の純冪零五代表。
全線形源はrootの固定源reportと
`CLASSIFICATION_DEGENERATE_QUADRIC_SOURCES_REVIEW_20261107.json` の独立全150再構成を入力とする。
以下で三源を除外し、残る二cofactor型を**三つの全既約源**へ帰着する。
§5の実変形によりnil31長鎖型IIは既知MPS成分の和に含まれる。
残る二源の成分帰属はまだ未解決であり、空と扱わない。

## 1. 直ちに除外できる三源

### nil31、radicalへの追加が零

q=2x0x2+x1²+x3²、X=x0∂x1−x1∂x2。
全十一源のEuler収縮Pfaffianは、Jacobiを課す前からx0qで割り切れる。
Poisson点ではこれが実退化式なので、rank四なら三次因子は可約である。

### nil31、短いblockからradicalへの追加

同じqでX=x0∂x1−x1∂x2+x3∂x4。
全十源のJacobiにはn5²と6n5n8−5n9の非零定数倍が現れる。
被約支持でn5=n9=0。この代入後、Euler収縮Pfaffianは全係数でx0qを因子に持つ。
従って既約三次因子の支持は空である。

### nil22、radicalへの追加が零

q=x0x3−x1x2、X=x0∂x1+x2∂x3。
U=(u,v,w)=(x0,x2,x4)上の三括弧は、全十五源で

    {u,v}=F_w, {u,w}=−F_v, {v,w}=F_u,
    F=(n8u²+2n12uv−n1v²)w/2
      +(n5u³+3n7u²v+3n10uv²−n0v³)/3

に一致する。全Jacobiには判別式n1n8+n12²と

    n0n8+n1n7+2n10n12,   n1n5−n10n8+2n12n7

が現れる。
`SMOOTH_QUADRIC_NIL_REVIEW.md` §3–4と同じCasimir論法が適用できる。
Nの全不変三次式はC3(U)+qℓ(U)であり、qからw²がなくなっても、この13次元空間は変わらない。
F=0では全五Pfaffianが零。rank四かつc既約ならℓ=0、cはFの非零定数倍でFも既約。
一方、上の三必要条件と同時SL2規格化により、binary quadratic係数がrank一ならu²|F、
rank零ならFはbinary cubicとなる。従ってこの全源も除外される。

## 2. nil31、長い鎖からradicalへの追加：二つの全二源が残る

q=2x0x2+x1²+x3²、X=x0∂x1−x1∂x2+x2∂x4。
全八源の変数をn0,…,n7とし、a=n2、b=n7、d=n6と書く。
Jacobiに現れるn5(3b+1)とn5(2b−1)からn5=0。
残る必要条件には

    E1=2ab+4a+2d²+2b²+7b+6=0,
    E2=4ab−2a+4d²+10b²+3b−4=0,
    d(3a+3b+2)=0,
    2n1d+3n4(2a+b+1)=0,
    3n0d+n3(9a+6b+5)=0

がある。E2−2E1を取ると

    a=(6b²−11b−16)/10,
    10d²+(2b−1)(3b+1)(b+2)=0,
    d(6b+7)(3b−4)=0.                                  (1)

これらの消去は多項式恒等式として検証した。

d=0の場合は(b,a)=(1/2,−2),(−1/3,−7/6),(−2,3)の三枝。
各枝でn3=n4=0、n0,n1自由で、全Jacobiは恒等零。
実退化式は順に0、(5/36)x3³q、(5/4)x3q²となり、全て対象から除外される。

d≠0の場合は次の二型だけが残る。n0,n1は任意である。

| 型 | a=n2 | b=n7 | d=n6 | n3 | n4 |
|---|---:|---:|---|---|---|
| I | 1/2 | −7/6 | ±5i/6 | −6dn0/5 | −4dn1/5 |
| II | −2 | 4/3 | ±5i/3 | 3dn0/5 | 2dn1/5 |

全てn5=0であり、各行の全二自由係数で全Jacobiが恒等零になることを確認した。
各±は射影反射x3→−x3とパラメータn1→−n1で交換される。
この反射の全bivector恒等式も確認したので、正符号を代表とする二つの全二パラメータ源でよい。

ここでは三次因子が実際に既約である。C=D/qと置くと、全自由係数について

    ∂C/∂x4=(d/2)x0²,
    [x1³]C=−d/6≠0.

Cはx4について一次で、その係数x0²と定数項の最大公約因子は1。
Gaussの補題からCは既約である。したがってこの二源を排除することはできない。

## 3. nil22、radicalへの追加非零：一つの全四源が残る

q=x0x3−x1x2、X=x0∂x1+x2∂x3+x1∂x4。
全十源のJacobiにはn9が現れるので、n9=0。
このときn0=0なら、Euler収縮Pfaffianは恒等的にx0qで割れる。
従って既約三次を要求する支持ではn0≠0。

Jacobiの条件n0(n5+n8)=0からn8=−n5。
さらに

    9n3n9−2n5²−n5n8+14n5+n8²−n8−20=0

は15n5−20=0となる。残る三つの必要係数で全てを解くと

    n9=0, n5=4/3, n8=−4/3,
    n7=−3n0n3/5,
    n4=−9n0n3²/25,
    n1=−27n0n3³/125.                                   (2)

自由係数はn0,n2,n3,n6、開条件はn0≠0である。
(2)を全括弧へ代入すると、全自由係数で全Jacobiが零になり、実退化式は

    D=−q [5x0q+n0(x2+3n3x0/5)³]/3.                     (3)

三次因子はx3について一次で、係数は−5x0²/3。
x0=0の定数項は−n0x2³/3なので、n0≠0なら係数との最大公約因子は1。
Gaussの補題により全開で既約である。
これは残余点の一例ではなく、同cofactor型の既約三次支持を尽くす全四パラメータ源である。

## 4. 有限証拠と次の義務

`code/classification_degenerate_quadric_q4_nil_support.py` はQQ(i)上のGaussian Fractionだけを使い、
producerとSymPyをimportしない。
全sourceの発散・q-normal・cofactor不変性、必要Jacobiのliteral係数、
全枝の普遍Jacobi、全実Pfaffianと五gradient、±反射を検査した。
`jobs/degenerate_quadric_q4_nil_support.json` のgate後、約0.93秒でPASS。
記録は `reports/CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_SUPPORT_20261107.json`。
全150線形源の階数は、前記の共通独立source reviewを依存とし重複計算していない。

全150係数での有限接空間計算
`CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_COMPONENTS_20261107.json` は、
I、II、nil22非零の順でJacobi微分の階数100、102、100、実族接の階数はいずれも47を与える。
従って実族の次元下界は22、射影接空間の次元上界は25、23、25である。
これだけでは成分性を主張しない。

IとIIの二自由係数は全て軌道方向である。元の複素座標で
κ=3(2a+b+1)と置くと、Iではκ=5/2、IIではκ=−5。
源の(n0,n1)=(0,0)の点を線形shear

    x4'=x4−n0*x0/κ−n1*x3/κ

でpushforwardすると、任意(n0,n1)の源と全係数で一致する。
上記reportはこの恒等式を全二自由係数のまま保存している。

nil31長鎖型には次のcofactor変形が用意できる。
(a,b,c,s,z)=(x0,x1,x2,x3,x4)と書くと

    q_t=2ac+b²+s²+5t²c²−4t⁴z²,
    X_t=(5t²b+4t⁴z)∂a+a∂b−b∂c+c∂z.

全tでX_t(q_t)=0。t≠0ではq_tはrank五、X_tの固有値は0,±t,±2tとなる。
これは**cofactorの変形だけ**であり、現時点ではΠ自体のPoisson変形を与えたとは主張しない。
既存正則源の座標引き戻しから全Πの極限を作るための入口である。

## 5. nil31長鎖型IIの全源はMPS閉包に含まれる

前節の軌道恒等式により、IIの原点パラメータの一代表だけを扱えばよい。
`x3_new=i*x3_old` と非零定数倍を施した有理代表Pは、
`CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_COMPONENTS_20261107.json` に全十括弧を保存した。
全150係数のkernelから実族接47方向を割った一つの真のnormal Rは

    R01=−6x0x2+x1²−4x1x3+x3²,
    R02=2x1x2−4x2x3,
    R03=−6x0x2−2x1²+2x1x3−2x3²,
    R04=−3x2², R12=3x2², R13=−2x1x2+4x2x3,

で、他の上三角括弧は零である。T01=x2²/5、T03=2x2²/5、他を零とすると

    P_t=P+tR+t²T                                          (4)

は全tで射影Poissonである。これは有限jetからの推論ではなく、
全350係数の射影Jacobi `E∧[P_t,P_t]=0` の多項式恒等式を直接検査した結果である。
strictなaffine Jacobiは零ではないので、P_tそのものをunimodular liftとは扱わない。

実退化式は非零定数倍を除きq_t c_tであり

    q_t=−t x2²+10x0x2+5x1²−5x3².

q_tは全tでrank四。c_0は§2の既約三次であるから、被約な既約三次の開性により
一般のtでもc_tは既約である。
divergenceをZ_tと書くと、標準unimodular lift
`Π_t=P_t−E∧Z_t/5` のcanonical q-cofactorは

    X_t=(2t x1/5−4t x3/5)∂x0
       +(6x0−8t x2/5)∂x1−6x1∂x2−4t x2∂x3/5+6x2∂x4.

`dq_t·Π_t=q_tX_t` を全係数で確認し、その特性多項式は

    λ(−5λ²+12t)(−5λ²+48t)/25.

t≠0では全五固有値が相異なり、二対の非零固有値の絶対比は2である。
従って一般のtはrank四の正則半単純cofactor型に属し、
`DEGENERATE_QUADRIC_CUBIC_REGULAR_SUPPORT.md` の全支持定理からMPS成分の和に含まれる。
その有限和は閉なので、t=0のPも含まれる。
最後に全二源が一つの軌道であることから、IIの全源の被覆が従う。

証拠は `CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_CURVE_20261107.json`
（全tの射影Jacobi、実五次式、0.34秒）と
`CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_CURVE_RECOGNITION_20261107.json`
（unimodular射影、cofactor、特性多項式、0.16秒）。
独立算術再構成は依頼中であり、成分の局所唯一性・正規性は主張しない。

残る義務はIとnil22非零の全支持の帰属である。
