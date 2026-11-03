# Vertex・scroll の有限交差と quad78有限欠損の独立監査

2026-11-03。`S61_NONSQUARE_VERTEX_SCROLL_FINITE_20261103_v63` の actual resultant、射影 coverage、全branchの有限性、corrected conditional への適用を独立監査する。旧omit-B²推論は使わない。有限G2の既存独立receiptはruntime復帰後も全hashが一致した。

## Scroll integrality の全点 coverage

元804項と四座標 frame は独立済みbinary reviewに固定する。actual二線形式の全8係数を読み、beta=`32a³J²/243`を確認する。第一行との比をzとすると、B=−zWで、(A,X,W)に対する二つのlinear equationを得る。betaがunitなのでrank低下は二つのquadratic r1,r2の同時零と同値である。

P²だけで分母を払い、固定4角Sylvesterを独自Fractionの24項determinantで全係数再計算する。結果は

`−67108864*a^9*J^6*P²*Eold/4749088069944112812890625`。

自然変数への多項式代入を直接行い、`Eold=P³E5`、

`E5=2(v+2)(7v−11)^2−t(55H−66v+18)^2`

を確認する。E5を全開の既知unitとはしない。

E5≠0なら全有限zで二行のrankは2。比を[s:tau]に斉次化すると、s=0でも二行の最初の二列のminorが40tau²/c²で非零となる。cross product(Ac,Xc,Wc)を用いた

`[s:tau]→[sAc:sXc:sWc:−tauWc]`

はbasepointのない三次morphismである。codeはuniversal linear coefficientsの多項式環でcross product、全minorへのfactorization、s=0でWc=40tau²/c²を確認する。

Cの第一行が非零な任意の点では比zが一意であり、rank2のkernel lineがその点を一意に復元する。第一行が零の点はB軸の一点だけであり、無限比の像である。従ってmorphismはCのsupportを尽くし、一般次数1となる。pullback O(1)=O(3)により像は次数3のirreducible curveである。一方Hilbert–Burch scheme Cも次数3、pureでCMである。従ってgeneric multiplicityは1、CMのS1とgenerically reducedからreducedとなり、Cは幾何学的integralである。

これはsmoothnessを仮定する議論ではない。別の同値な確認として、第一行非零のopenのschemeはrank2 kernel sectionによりA1と認識でき、CM/S1がB軸点でのnilpotent/embedded成分を排除する。

## Vertex比例とE5退化の有限交差

vertex binary matrixの列は正確に[2,4,5]。選ぶ二つのminorは[2,5]と[4,5]であり、元sourceからの正規化はQQ定数とtのunitだけである。B*S/B*Uを比較する旧推論は関係しない。

E5=0をD²t=Nと書く。D≠0ではt=N/D²が一意であり、二つのlinear-t必要条件からF_i=D²beta_i+Nalpha_iを得る。保存された正規化はQQ定数だけで、D,J,N,fの因子を除いていない。全係数を独自Fractionで確認する。

F/GのH係数に対するQQ[v]和1と、v=0でH次数3,2を保つQQ[H]和1を確認する。共通vertical factorは前者で排除される。H正次数の共通factorは両leading coefficientを割るので、この特殊化でも正次数を保ち、後者に反する。よって共同零点は有限となる。

D=0ならN=0であり、因数式から(v,H)は(-2,−30/11)または(11/7,120/77)の二つだけである。保存リストがこの二枝を正確に尽くすことを確認し、actual二つのlinear-t式と保存QQ[t] cofactorの和1をそれぞれ独立検査する。両枝は空である。従ってE5=0とvertex rank低下の交わりは全aPJ-openで零次元以下。

## G1の有限coverは全aPJ-openで使える

既存proportionality独立receipt自体がaPJ-openでの排除を証明している。DB≠0の係数unitは−a³J/12。DB=0で用いるL=20P+3e+20は−180a/Pなのでunitであり、cross-minorの除算もP,L,Jだけである。F,B0,Cの非零性は使わない。

G1 failureが残るにはQ0 mod(B,W)がrank2でなければならず、rank4のsplit formからDB=0、dual norm N=0が必要となる。共通linear factorを持つQ3もrank≤2なので、その保存3角minor M=0が必要である。座標とminorの正規化はa,P,L,Jだけを使う。ここもF,B0,Cを仮定しない。

既存のfinite-Sigma独立監査はNのe⁴係数が非零定数、P²⁸多項式が(N,M)に属し、そのleading coefficientも非零定数であることを確認した。この二つをmonic化すればQQ[P,e]/(N,M)自体が有限次元になる。後段のFULLU-unit・有限体・Nakayamaによるraw24全射性は、この有限性には不要である。従ってG1 failureは**全aPJ-open**で有限coverに入る。

## Corrected criterion と次数3 colon

G2 failureは先行の独立ternary-pencil有限性で有限。E5≠0ではCがintegral、E5=0でも今回の有限交差を除けばvertex rank2なので、corrected conditionalのG1∧G2∧(G3またはintegral C)が適用できる。従ってraw84×126のrank78欠損は全aPJ-openで零次元以下となる。有限集合上での実rankは未検査である。

元light rank146への移行ではfull CIを仮定しない。四式のidealをJ0とし、次数4の関係ell*f3=ΣQi*Ai2をmod ellへ落とす。raw rank78によりAi2 mod ellは六Koszulの定数結合のみなので、そのactual liftsを引いて全Ai2をell倍にできる。ellを消せばf3∈(J0)_3、すなわち `(J0:ell)_3=(J0)_3`。codeはこのKoszul-lift cancellationをuniversal polynomial identityで確認する。

六式I=J0+ell(Z,W)について必要な次数で `(I:ell)_3=(J0+(Z,W))_3` が得られる。raw degree4商は126−78=48、元ZW商のdegree3は35−19=16なので元degree4商は64、light像rankは210−64=146となる。六式の独立性とheavy列のdisjoint性からweighted rank152も従う。これらの元chartへの適用はFULLU上、有限集合を除いた範囲である。

有限集合の空性、全点rank78/lower146、full cubic57、残る全quartic quotientの全射、長65、ambient除去は未証明である。誤ったomit-B² syzygy適用は復活させない。

再現（新しい出力先）:

```sh
source /workspace/cloud-setup/activate.sh
python code/s61_nonsquare_vertex_scroll_finite_review_20261103_v63.py \
  --output work/marked-square-v63/vertex-scroll-finite-independent-new
```

60秒・2GB以内の係数恒等式照合だけを行い、GB、新しい体、prime、重い行列探索は起動しない。
