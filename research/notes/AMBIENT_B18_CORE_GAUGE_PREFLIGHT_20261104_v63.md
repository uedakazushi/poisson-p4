# Tiny core の実group cutsとalgebraic rank2 recipe

2026-11-04。相対weight≤3 theoremのfrozen bytesは変更しない。同じ実CA/CB core10（105/127項）と同じG19切断だけを調べる。新しい点、source再生成、GB、generic coefficient solverは使わない。

## 1. Coreはまだ共通G19 sliceにない

十九切断のcore上の中心値からの差が非零なのは、両者ともraw rows

    15,16,22,23,50,65,132

である。全polynomialをcertificateに保存する。row50はCAで−2*k2*p、CBでk2*p/2なのでfixed零weight sector(k2=0)では六本になる。

共通zero baseにしたい九raw rowsは

    27,108,109,112,113,117,118,119,61。

**group gaugeを掛ける前の**九函数は、p以外の八zero source変数K8に関してlinearだが、p=0でのrankは7しかない。従って、この九函数をそのままcommon R0と扱い八linear pivotを逆にする方法は使えない。QQ(p)で非零になるminorでもpのfactorはfold parameterなので単元にできない。保存した九minorとgcdはこの点を明示する。

## 2. Constant Schur補正はunit八pivotとfoldを回復する

固定torusに関して十九group方向中のzero15と、zero切断15を取る。対応する定数matrix Qは可逆である。保存O19とcoreのK8微分Lを用い、linear gaugeは

    gamma_linear=−Q^-1 L_cuts K,
    L_base_corrected=L_R9−O_R9 Q^-1 L_cuts

と正確に求まる。後者のrankは8であり、八raw base rowsとunit minorを保存する。残るrowがfold方程式となる。これは高次gaugeをlinear式で置き換えたものではなく、次のimplicit systemのJacobian selectorを決める計算である。

同じsource coreの既存p²係数WをこのSchur操作で投影する。さらに八base方程式を解いた後の一つのfold quadratic係数も非零となる。全23×23 blockのdeterminant、八pivot、残るfold rowと正確な二次係数を保存する。

## 3. 未展開だが有限なpolynomial implicit recipe

前packetと同じgroup sectionを

    G=I+sum_{j=1..18} gamma_j E_j, tensor multiplier=1+gamma0

とする。quadratic bivector Uのpushforwardの分母はdet(G)²であり、分子は

    (1+gamma0) G U(adj(G) X) G^T。

選択raw row cごとの十九gauge equationは、この分子のc係数からP0_c det(G)²を引く。これは明示的な有限polynomial circuitで、group変数次数は高々11。巨大展開しない。det(G)、1+gamma0、gauge Jacobianは中心で単元として扱う。

zero sectorではsource k2=0とpositive group4を零にする。equivarianceとimplicit解の一意性により、positive group4は全次数で零である。残る15gauge equationsに八選択base方程式を加えると、未知gamma15,K8に関する23次Jacobianが単元となる。したがって九target変数xとpの上に、中心gamma=K=0を通る唯一の**algebraic power-series** solutionが存在する。ここでformal IFTをfinite polynomial systemへ適用しているためalgebraicityが得られ、単に微分を切ったものではない。

残る一つのbase equationをF(x,p)=0と書く。中心でp一次は零、保存した二次係数は単元なので、pointed Weierstrass factorはdegree2である。従来の全source finitefree2 theoremとも一致し、完成normalizerは

    QQ[[x9]][p]/(p²+a(x)p+b(x))

となる。a,bの実係数は未展開である。

algebraicityも有限recipeから得られる。23unit equationsと最後のfoldは局所的にquasifiniteなalgebraic germを定めるので、その座標はQQ(x9)上algebraicである。二つのlocal formal conjugatesのsum/productがa,bとなり、algebraicかつformal baseに属する。各source座標のbasis(1,p)係数も二共役のtrace/differenceからalgebraicとなる。従ってtracefree quadratic normalizerとactual target係数は、algebraic power-series係数で構成する有限rank2 recipeを持つ。

これは**完成環またはalgebraic-power-series環上でのfinitefree2**であり、普通のpolynomial base上のfinite/proper modelを得たとの主張ではない。unit-Jacobian algebraic implicit branchだけから未完備のno-escapeを推定しない。rational/triangularなclosed GL gaugeを発見したとも言わない。

次の実計算対象は、このselectorで指定したalgebraic gaugeとquadratic Hensel factorの係数、次いで共通zero17 targetへの二つの実像と六cubic graph差である。既存full Dのnil/socle/depthはまだ決めていない。
