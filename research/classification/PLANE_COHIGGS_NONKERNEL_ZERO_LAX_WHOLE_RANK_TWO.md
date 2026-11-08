# Co-Higgs M(p)非零・T=R零の全四源：定数rank束でrank二以下

2026-11-08、PRIMARY E、log_manuscript。producer純幾何証明。
新しい計算、M1、minor、point、profile、GB又はjetは実行しない。
独立受理は別義務。1054主稿／frontier／global／Gitと旧bytesは不変。

## 1. 正確な四whole source

rank零planeのO(1) co-Higgs、rank M二、M0=M(U)非零、
T=R=0の全sourceを扱う。
[旧全high-m入口](PLANE_RANK_ZERO_RANK_TWO_HIGH_MULTIPLICITY_FIRST_JET.md)による全actual normalizerを使うと、
これらは正確にNDSS0、NDN0、BSS0(tau=0)、BN0(tau=0)の四つである。
残るBSS0／BN0のtau非零は本命題に入れない。

M:C3→sl2、n=(z,w)、q_A=det(n,An)とすれば、全first symbolは

~~~
P0=U∧(M(x)n·∂n) − q_(M(U))(n)∂z∧∂w/2 .
~~~

これはT=R=0、S=−M(U)/2の全Lax解であり、全multiplicityで有効。
全P1三十及びFULL九P2を保つ。

## 2. 全first-symbolの実rank三への変形

H=diag(1,−1)、E=E12、F=E21とする。
各元sourceでUとM0を正確に保つ次の一本のsymbol familyを取る。

|target whole source|固定U|M_epsilon|固定M0|K|
|---|---|---|---|---|
|NDSS0|∂b+∂c|epsilon aH+bF+cE|F+E|−(z²−w²)/2|
|NDN0|∂b|epsilon aH+bF+cE|F|−z²/2|
|BSS0(0)|∂b|epsilon aF+bH+cE|H|zw|
|BN0(0)|∂c|epsilon aF+bH+cE|E|w²/2|

T=R=0とするので、各epsilonで全Lax式は恒等的に成立する。
epsilon非零では三つのcoefficient matrixが独立、rank M三。
M0は全parameterで非零かつ同じss又はnil型である。
従って[rank M三tau零の全源定理](PLANE_COHIGGS_FULLRANK_HIGH_M_WHOLE_RANK_TWO.md)が
epsilon非零の全higherへ適用できる。
epsilon=0は元の四whole sourceそのものである。
元P1だけを固定してこのsymbol curveへ代入するのではなく、次節で全higherを正則に延長する。

## 3. 元の全kernelを束として延ばす

D1(epsilon)は元全350×30、D2(epsilon)は350×9の写像で、
P0(epsilon)の係数がpolynomialなのでentryもpolynomialである。
以下はその新しいoperatorを計算する手順ではなく、既存rank定理の幾何的適用である。

epsilon非零では前掲全rank三定理によりrankD1=22、rankD2=9。
epsilon零の四targetについては旧全源独立reviewが既に次を確認している。

|target|rankD1／全kernel|既存22minor|rankD2／全kernel|既存9minor|
|---|---|---|---|---|
|NDSS0|22／8|2187/512|9／0|6|
|NDN0|22／8|6561/4096|9／0|9/2|
|BSS0(0)|22／8|11664|9／0|−288|
|BN0(0)|22／8|−19683/512|9／0|−27/2|

元report CLASSIFICATION_PLANE_RANK_ZERO_RANK_TWO_HIGH_M_20261107.jsonの
NDSS0_reused.old_full_certificate、NDN0.analysis、
BSS0_tau／BN0_tauのparameter零certificateをそのまま読む。
新しいrank／minorを計算したのではない。
全値と全350行の独立認証は
[all-m独立review](PLANE_HIGH_M_RANK_AT_LEAST_TWO_SEPARATE_REVIEW.md)を入力とする。

旧22minorと旧9minorが非零となるepsilon零を含むopen Bへ縮める。
epsilon非零で全23minorが零なので、多項式恒等式として全Bで零。
旧22minor非零と合わせ、D1は全Bで定数rank22、
K=kerD1はrank八vector bundleである。
D2も全Bで単射で、旧九rowに対するinverseはB上正則である。
よって

~~~
Q_epsilon(s)=−(D2_selected(epsilon))^(-1) Jac(s)_selected
~~~

は全s∈K上の正則な二次mapとなる。
任意の元full Poisson tensorは必ずP2=Q_0(P1)を満たす。
選択九式で満たさない他の元式は捨てず、次節で束全体上の恒等式として戻す。

## 4. 一般fiberの全sourceが本当に稠密である

rank三tau零の前掲証明では、全八kernelは有限shear六＋二方向N1,N2へ移り、
任意の二係数についてJac(N1,N2)=0、全P2=0である。
従って一般fiberの**全K1**がactual Poisson tensorへ延び、その全actual Pfaffianが零となる。
「そのfiberのPoisson subsetだけがrank二」という弱い結論を使っていない。

全Kはintegral base B上のvector bundleなので既約、epsilon非零の全Kは稠密である。
P0(epsilon)+s+Q_epsilon(s)の全350Jacobi係数と全126actual Pfaffian係数はK上正則函数。
一般fiber全体で零なので、全K上で恒等零。
特にepsilon零の全八higher係数でもPfaffianは零である。

これは任意の元Poisson solutionが入る全候補を覆う。
K0のchosen slice、一般higher点又は基底だけの閉包ではない。
定数rank八の束に全P1を保持し、FULL九P2の唯一正則復元を使うため、
rank-jumpから新しいhigherが生じる補集合を残していない。

従って四whole sourceの全projective tensorはrank二以下。
元のactual normalizer逆、全zero係数、全higherを戻して結論する。
全BSS0(tau非零)、BN0(tau非零)や他co-Higgs sourceの空性は本稿では未主張。

## 5. 証明の射程

本証明は新しいCAS又は算術reportを必要としない。
一般fiber全sourceのrank二定理、既受理の四targetの完全kernel／D2単射、
explicit M_epsilonとvector bundleの正則性だけからなる。
これらの独立読解を受ける前にglobal claim又は主稿へ採録しない。
