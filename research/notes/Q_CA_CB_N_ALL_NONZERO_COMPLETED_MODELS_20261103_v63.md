# 全非零Qaで同じCA/CB/N完成環：Q1だけのpinchではない

2026年11月3日。対象は個々の被約成分CA,CB,Nの局所環である。
全Poisson局所環の共鳴a=2/9や全枝の分類とは区別する。

## 結果

全ての非零複素aについて、物理的なQaにおける被約CA,CB,Nの完成局所環はそれぞれ

    CAhat ≅ C[[u1,...,u19,s1,...,s8,p²,p³,s1*p,...,s8*p]],
    CBhat ≅ C[[u1,...,u19,s1,...,s8,p²,p³,s1*p,...,s8*p]],
    Nhat  ≅ C[[u1,...,u19,s1,...,s7,p²,p³,s1*p,...,s7*p]]

である。

|target|dim|embedding dimension|depth|multiplicity|
|---|---:|---:|---:|---:|
|CA,CB|28|37|20|2|
|N|27|35|20|2|

従来のa=1の型は、同じ非零Qa曲線全体に適用される。
a=2/9でもこれら三つの被約成分の完成環は変わらない。
全Poissonスキームにおける共鳴方向の増加は、これら個々の成分の局所型の変化を意味しない。

a=0での物理的Qaはrank2であり、上の完成型を延長していない。
ただし後述の閉性から、CA,CB,Nは物理的Qa0でもnonnormalかつnonCMである。

## 1. 正確なS61曲線と、動いているのはfree k0であること

既存の実物理normalizationと重み付きscalar作用により、a≠0のQaはS61内の

    u25=−4/3, u36=3, u38=9a/2,
    その他の全61座標=0

に移る。scalar normalizationのframeはaを逆にするため、物理的a=0にはこの同定を使わない。
このS61内の式自体はa=0にも延長するが、その点は元の物理的rank2 Qa0ではない。

accepted closed proper sourceの座標で、この全曲線は

    beta=(tau,A,B,C,D,F)=(0,0,4/3,3/2,0,0),
    r=0, k0=9a/2, k1=...=k8=0

である。fiber parameter順序は(r,k0,...,k8)であり、k0とk1を混同しない。
βは7次元critical rank2軌道上で固定され、変化するk0は既存pinched-imageの**自由な三fiber方向の一つ**である。
欠損odd idealが(a−1)になるという仮説は成立しない。

全S61曲線で原134式、150係数再構成、全projective Jacobiを厳密に確認した。
CA,CBの閉source式への同じ代入も全aの恒等式として一致する。

## 2. 全aでの実際の旗fiber

reference lineを含む旗chartで

    g=[[1,r,−w,z],[0,1,0,0],[0,z,1,0],[0,w,0,1]]

を使い、全61座標を有限に変換してclosed CA/CB constraintsへ代入した。
その全多項式をQQ[a,r,z,w]で計算した結果、CA,CB,Nの全旗idealは正確に

    (r,z,w²)

である。QQ(a)でのgeneric計算だけでなく、aを多項式変数に含めたGroebner basisで同じidealを得た。
従ってaに関する例外分母・未検査の特殊値はない。
K_betaは固定で、そのkernelはspan(e_y,e_w)である。
chart外の唯一の候補e_wはaに依存しない基底constraintで除かれるので、これは**全projective fiber**である。

全fiberは一点の台を持つ長さ2である。
proper smooth sourceはその近傍で有限normalizationとなり、targetに同型ではない。
この事実は全非零aでnonnormalityを直接与える。

## 3. 全pinched-image theoremの適用

既存の全像定理 `AMBIENT_COMPONENT_PINCHED_MODELS_20261103_v63` は、
同じcritical β軌道上で、actual target q座標Kを使い

    A=R+I*p+p²*Bnorm

を全次数で証明している。missing idealは

    N:  I=(s,k1,k3,k4,k6,k7,k8),
    CA,CB: I=(s,r,k1,k3,k4,k6,k7,k8).

ここでk0,k2,k5は欠損idealに入らない。
基底をB=C=1へ戻す固定Sp4 torusはmonomialをdiagonalにrescaleするだけなので、
k0だけが非零という全Qaの形を保つ。
既存証明書に保存された全45 odd係数へ、k0を変数としてその他を0とする代入を行い、
三族全てで厳密に0であることも確認した。

従って全Qaは同じnonnormal stratum内にあり、既存の全像環Ac8/Ac7が全ての非零aに適用される。
九つの値・perspective smooth方向を戻して上の全射影完成環を得る。
これはsource微分の階数だけから完成環を推測する議論ではない。

normalizer内のconductorは(s1,...,sc,p²)、target内では
(s1,...,sc,p²,p³,s1*p,...,sc*p)である。
normalization quotientはC[[u19]]*pで、depth20などの数値はこの正確な列から従う。

## 4. 物理的a=0で残ること

独立した全QQ[a] Smith証明書 `Q_CA_CB_N_SOURCE_SMITH_20261103_v63` は、
物理的proper-source像微分のaffine rankが

    CA: 28 全a,
    CB: 28 (a≠0), 27 (a=0),
    N: 27 全a

であることを証明した。これらはtarget tangent dimensionではない。
特にCBのa=0のsource geometryは別に調べる必要がある。

excellent reduced varietyのnormal locusは開である。
全a≠0のQaが三成分のnonnormal locusに入るので、その閉包のQa0も各nonnormal locusに入る。
三成分それぞれの既存の**全射影R1定理**により、Qa0での局所環もCMではない。
局所環としてのSerre S2が成り立てばR1と合わせてnormalになるためである。
これはclosed pointのdepthが2未満という主張ではない。
物理的a=0でのfiber次元・normalization ring・全完成型はこの閉性だけでは決まらず、未解決として残す。

## 証明書

* `code/q_CA_CB_N_all_nonzero_completed_models_20261103_v63.py`
* `reports/Q_CA_CB_N_ALL_NONZERO_COMPLETED_MODELS_20261103_v63/report.json`
* `code/q_CA_CB_N_source_smith_20261103_v63.py`
* `reports/Q_CA_CB_N_SOURCE_SMITH_20261103_v63/{report.json,certificate.json.gz}`

新規計算は全a flag/image検査2.87秒、全physical source Smith1.84秒。
全Poissonのfull35完成環・残る枝・実CA/CB scheme交叉の全Qaでの厚みは、この結果から分類していない。
