# rank四の二次因子：nil(3,1)・nil(2,2)五源の全支持還元

2026-11-07。課題 `E_reduced_nonnormal`。
対象は `DEGENERATE_QUADRIC_CUBIC_COFACTOR_REDUCTION.md` のrank四の純冪零五代表。
全線形源はrootの固定源reportと
`CLASSIFICATION_DEGENERATE_QUADRIC_SOURCES_REVIEW_20261107.json` の独立全150再構成を入力とする。
以下で三源を除外し、残る二cofactor型を**三つの全既約源**へ帰着する。
§5の実変形によりnil31長鎖型IIは既知MPS成分の和に含まれる。
§6の真三方向変形と忠実graph認識によりnil31長鎖型IはT40に含まれる。
§8の滑らかな二次曲面への実変形によりnil22非零の全源もMPS成分の和に含まれる。
従って対象の五cofactor型は全て除外又は既知成分の和に被覆される。

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
`CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL_CURVE_REVIEW_20261107.json` は
producer/SymPyを使わない独立Fraction実装で、全tの射影Jacobi、unimodular gauge、
全五cofactor式、120項の特性行列式、実五次式の因数分解を再構成した（0.014秒、PASS）。
成分の局所唯一性・正規性はこのIIの証明には不要であり、主張しない。

## 6. nil31長鎖型Iの全源はT40に含まれる

Iの有理代表Pは座標重み(2,1,0,1,−1)、bivector重み−1で同次である。
全150係数でgroup・Euler-vertical・scalar接方向47次元を除いた安定線形sliceでは、
真のnormalは三つで、その相対重みは4、2、3となる。
全ambient相対重みの範囲は−4から6なので、同変陰関数の重み付き次数上界は6。
有限計算は全補正項が零、全残余Jacobi idealが零となることを示す。
従ってliteral三方向R0,R1,R2に対して

    P+aR0+bR1+cR2

は任意(a,b,c)で実際の射影Poissonテンソルである。
`CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL31_I_WEIGHTED_20261107.json` に
三方向・安定slice・全100独立Jacobi行・全族を保存した。
この族の三方向は全kernelをgroup接47方向と合わせて張るため、
実族の射影次元下界は25、全接次元上界も25である。
ゆえにPは全射影Poisson概形のsmooth25点で、一つの成分だけを通る。

同じPを点(1,0,1,0,1)でDarboux規格化して忠実S61 chartへ移し、
全150係数の商再構成・全134方程式・T39/T40の全45graph式を検査した。
結果はT39には属さず、T40の全被約graphに属する。
記録は `CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL31_I_GRAPH_20261107.json`。
したがってPを通る唯一の成分はT40である。
§4の全二源の軌道恒等式から、Iの全源がT40に含まれる。

`CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL31_I_REVIEW_20261107.json` は
独立Fraction実装で全D100、group47、三実方向を加えた階数50、全350の普遍Jacobi、
同一源点、全PGL変換と150係数の商同値、全134方程式及びT40の全45graph式を再構成した
（3.396秒、PASS）。
なお正規五次式を使うために行った二つのMacaulay rank診断は十分なrankを与えず、
この帰属証明の依存には一切使っていない。

## 7. nil22非零の全源も一軌道へ還元できる

§3の全源でn0≠0とする。λ≠0をn0=λ^(−2)となるように選び、k=3n3/5と置く。
原点代表(n0,n2,n3,n6)=(1,0,0,0)から、次の線形変換を順に施す。

    G0=diag(1,1,λ,λ,1),
    Gz: x4'=x4+3(n2−kn6)x0/10+3n6x2/5,
    G3: x2'=x2−kx0, x3'=x3−kx1.

合成G3 Gz G0によるpushforwardは、全四自由係数を保持した§3の源に一致する。
`CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL22_CANONICAL_20261107.json` はこの全十括弧の恒等式を検査した。
標準代表は重み(2,0,1,−1,−2)、bivector重み−2で同次である。

この代表の全150係数での有限同変sliceは、真normal三つの相対重み1、2、3を持つ。
全陰関数の次数上界は9であり、その全多項式と残余idealは
`CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL22_WEIGHTED_20261107.json` に保存した。
normal変数を(a,b,c)とすると、残余idealは非零定数倍を除き

    (11a³+18ab−54c) · (4a²+9b, ab+6c, 3b²−8ac).

従って被約sliceの支持は曲面11a³+18ab−54c=0と
曲線b=−4a²/9、c=2a³/27の和である。
これは任意精度jetからの外挿ではなく、全同変陰関数の有限上界内での全残余idealである。
ただし各枝の既知成分への帰属はこの計算だけでは従わない。
全源は一軌道なので、原点代表への一つの実変形を既知成分へ認識すれば全被覆できる。
この支持分解は別担当も独立読解して受理した。
Hは素ideal、Fはcに一次で滑らか、H上のFは−a³なので二枝は相異なり原点だけで交わる。

## 8. nil22非零の全源は滑らかな二次曲面への実変形で被覆される

全有限sliceの多項式に

    (a,b,c)=(t,−4t²/9,2t³/27)

を代入して得たP_tの全十括弧を
`CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL22_CURVE_H_20261107.json` に保存した。
このP_tは全tで射影Poissonである。
実際、ここでも有限jetではなく、literal多項式族について全350の射影Jacobi係数が恒等零であることを直接検査した。
その実退化式は非零定数倍を除きq_t c_tであり、二次因子は

    q_t=−24300x0x3+24300x1x2+4050t x1²+675t² x1x3
        −180t³ x1x4+72t³ x3²+2t⁵ x4².

五変数Hessianの行列式は

    det Hess(q_t)=1394713760400000000 t⁵.

従ってt≠0でq_tは滑らかである。t=0の三次因子は§3のn0=1の既約三次なので、
一般tでもc_tは既約である。
滑らかな二次曲面＋既約三次曲面の全支持定理から、一般tのP_tはMPS成分の和に含まれる。
有限個の成分の和の閉性によりP_0も含まれ、§7の全四自由係数の軌道恒等式から
nil22非零の全源が被覆される。

この最終被覆には有限sliceの全支持性や各枝の局所次元は不要であり、
全源の軌道恒等式、保存したliteral実曲線、実退化式の因数分解、
Hessian行列式、及び既存smooth-quadric全支持定理だけを用いる。
`CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL22_REVIEW_20261107.json` は独立Fraction実装で
全十括弧のG3GzG0 Laurent普遍恒等式、det G=λ²、n0=λ^(−2)による全開被覆、
curve endpointの一致、c0のprimitive既約性、全tの射影Jacobi・実因数分解・Hessianを再構成した
（0.034秒、PASS）。
さらに `CLASSIFICATION_DEGENERATE_QUADRIC_Q4_NIL22_CURVE_REVIEW_20261107.json` も
別の独立Fraction実装でliteral curve、実因数分解、Hessianの全120項を検証した（0.022秒、PASS）。
