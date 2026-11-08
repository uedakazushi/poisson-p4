# Cubic-cone rank二の残critical源：T38と一つのactual24軌道閉包へ

2026-11-08、PRIMARY E。保存済み全kernelと新しい全四critical係数報告からのwhole-source純証明。
A_ss/A_nil各λ=2/5,3/5の全六fiberについて、
**λ=2/5はT38、λ=3/5は一つのactual二十四次元PGL軌道の閉包へ全点が入る**と示す。
後者の全ambient componenthood、既知親との同定又は新componentは未証明。
有限報告のsymbolic factor=[5]だけから絶対既約性を推論しない。

## 1. 元四全sourceの有限 coefficient 報告

`CLASSIFICATION_CUBIC_CONE_RANK_TWO_FOUR_CRITICAL_SOURCES_20261108.json` は、
旧complete27actualbasisと保存済み四exactwallの全6kernelを直接戻した。
全27source coefficient／十actual bracket／全350Jacobi／canonicaldiv零／全126actualPfを確認し、
全actualDは旧sourceのparticular+vertexDmatrixと一致した。
activegate後2.77841秒、RSS51,792KiB。新M1、operator、λ生成、rank、minor又はnativeは無い。
四全sourceのpolynomialfactorは[5]で、properfactorによるwhole空性は得ていない。
以下はその全式からの新しい幾何的含意である。

q0=x0x2−x1²、q1=x0x3−x1x2、q2=x1x3−x2²。
座標cは保存報告のx4である。

## 2. λ=2/5の全semisimple六fiberは既rank零T38 sourceに等しい

旧rank零L_semのrawbaseは P02=q0,P03=2q1,P13=q2、他零。
そのcolumn-divergenceは

`d0=−3x0, d1=−x1, d2=x2, d3=3x3`。

五変数全unimodularcanonical liftは rawΠ−(1/5)E5∧div(rawΠ)。
元全六括弧を比較すると、A_ss(2/5)のcanonicalbaseは
**このrawL_sem baseのcanonical liftの−一倍**である。
例えばP01=2x0x1/5、P02=−x0x2/5+x1²、P03=−4x0x3/5+2x1x2で、
残三entryも同じliteral identityを満たす。
c-divergence項も同時に一致する。

旧rank零の全六vertex rawsource W をcanonicalへ写すlinear mapは

`W -> W−E*div(W)/5`。

これはinjective。kernelならW=E*f、fはlinear。
各W_iは(q0,q1,q2)内なので、twisted cubic cone上で全xi*f=0。
このintegral coneはP3をspanするためf=0、W=0。
旧wholekernelは六次元、新A_ss(2/5)のfullvertexkernelも正確に六次元。
canonical mapは全Poisson／全conormal sourceを保ち、baseも上の−一倍で一致するから、
これらのfixed six-dimensional vertex spacesは等しい。
固定したconstant linear同型の逆で全coefficientを戻せ、parameter unitの隠れた補集合は無い。

したがって全A_ss(2/5) sourceはEuler／nonzero scalarを除いて旧rank零L_sem sourceそのもので、
既 `CUBIC_CONE_ZERO_T38_WHOLE_RECOGNITION.md` と独立receiptからT38へ入る。

## 3. 同じcriticalλでnil全fiberへ戻すconstant-rank橋

各λ=2/5又は3/5を固定する。
旧fullbase支持のAbranchではnonzero binary quadratic Uに対しW=λU、Vbase=0。
semisimpleとnilpotent両actualSL2orbitでfull350×18vertexmatrixのrankは十二、全kernelは六。
したがってirreduciblebase C3\{0} 上で全kernelはrank六vectorbundle。
semisimple逆像はopen denseで、全totalspaceはirreducible。
全actualtensorへのmapはregularである。

よってλ=2/5のnil全六fiberもclosed T38へ入る。
λ=3/5のnil全六fiberも、以下で定めるsemisimple全sourceのPGL飽和閉包へ入る。
genericbasisのdenominator代入や未証明flatnessを用いず、全rank theoremに依存したvectorbundleの議論である。

## 4. λ=3/5：全kernelから二cohom係数への実shear逆

全sourceactualDは、保存六kernel座標 k0,…,k5に対し正確に

```
F=(c+ell)*q0*q2 + a*q2*Cminus + b*q0*Cplus,
ell=(5k4/2)x0−5k1*x1−(5k2/4)x2−(5k3/2)x3,
a=k0, b=k5/2,
Cminus=x2*q2−x3*q1,
Cplus=x0*q1−x1*q0.
```

四つのq0q2*xi方向は独立。
残二方向もそれらと独立である：q0で割った余りからa=0、q2で割った余りからb=0が従う。
具体的に

`x0*Cminus=−q1²−x2²*q0, x3*Cplus=q1²+x1²*q2`

であり、q0,q1,q2は相異なるirreduciblequadrics、各linearxiともcoprime。
従ってfixed全六vertexkernel→actualDのlinear mapはinjective。

actual c'=c+ell(x) はbasePをliteralに保ち、fullvertexkernelを保つ。
actualDにおいてellを消し、a,bを変えない。
上のinjectivityから、変換後tensorはk0=a,k5=2b、他四係数零のcanonical tensorへ戻る。
全coefficient上のactual inverseであり、Dだけを合わせてtensor同値を仮定したのではない。

## 5. 全unit／zero補集合と一軌道正常化

a=0ならFはq0、b=0ならq2を因子に持つ。
これらの全closed sourceは対象既約quinticでは空。
ab非零なら、cのcoefficient q0q2と残項のgcdは一。
上の二identityからq0はCminusを割らず、q2はCplusを割らない。
よってprimitive c-linear polynomialとして C上絶対既約である。

actualSL2Cartanのbase map diag(t³,t,t^{-1},t^{-3})は固定Pbase(3/5)を保つ。
vertex相対scale c→gamma*cと組み合わせると

`(a,b) -> (gamma*t^5*a, gamma*t^{-5}*b)`。

t^{10}=b/aとgamma=(a*t^5)^{-1}を選べば両方一。
全ab非零でcomplexrootを選べ、逆は非零diagonalの逆。
全rawcanonical sourceは固定baseとdivergence-zeroを保ち、D map injectiveなので、
このD係数の変化はactualtensorの同一正常化そのもの。
したがって対象全sourceは

`λ=3/5, k0=1, k5=2, k1=k2=k3=k4=0`

という **一actual PGL軌道**へ戻る。
新しい任意point/profileを選ぶのではなく、whole-sourceの全unit逆から得た代表である。

## 6. 実軌道dimension二十四

Fは(q0,q1,q2)²内なので、irreducible cubic-cone SはSingFに含まれる。
Fは被約既約quinticだからSingFはdimension高々二、Sはそのirreduciblecomponentである。
connected tensorstabilizerは有限個のsingularcomponentsをそれぞれ保ち、従ってSを保つ。

固定λ、semisimplebaseのdimensionは三、全vertexkernelは六で、fixedSaffinesourceはdimension九。
canonical27sourceはEuler商と全係数上で一対一なので、この九はactualtensor familyのdimensionである。
Aut(S)八次元とsource scalar一次元の群が、全baseSL2／scalar正常化、
四vertexshear、残二toruscharacterにより上のab非零openへtransitiveに作用する。
群もdimension九なのでstabilizerは零次元。
connected projective tensorstabilizerがSを保つことから全PGL5 stabilizerも零次元。
よって **actual PGL5 orbit dimensionは24**。

全semisimple sourceのab非零openはdenseで、その閉包は一つの24orbitclosure C。
a0/b0の全補集合もその閉包に入り、§3からnil全六fiberもCへ入る。
これは全sourceのclosed-family認識義務を一つの実軌道閉包へ縮めるが、
Cが全ambientcomponentか、高次元既知componentのboundaryかはまだ分からない。

次の最弱義務はこの同一代表一つのfullambient D／Euler／実群作用を確認すること。
新しい点の列、rankprofile、jet又は未知parameterのsampleへ進まない。
