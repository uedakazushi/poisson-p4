# CA実self-pair族のrank2・二つのimmersion・cross-tangent

2026-11-03。独立監査済みの実QQ[a] self-pairを入力とし、接空間計算と全proper fiberから従う局所結論を区別する。
対象は個別のreduced CAであり、full Poisson schemeやCA/CBのactual intersection Dではない。

## 実テンソルとsourceの微分

先行packetのP(a)=∂y∧V(a)は全aでPfaffianが零である。
(Euler∧P)012のx³係数は1なので、射影Poisson rankは**全aで正確に2**である。
generic-rank4のself-pairは構成していない。

両旗liftをfixed flagへ運んだテンソルをPとQ=g_*Pとする。
実CA fixed20行列Cと全flag9の物理係数Lie微分から、両sourceの全150×29微分を再構成した。
Cの20次単元blockを消去すると、二つの150×9行列のrankは9であり、選択minorは

    first:  −11943936*a^9/15625,
    second: 864*a^14/15625.

従って両source微分は**全a≠0でrank29**である。
単一標本からの外挿ではなく、全QQ[a]行列とQQ(a) row reduction、厳密なminorである。

## Orbitと同じ点における二接空間

実25個のGL Lie微分のQQ(a) rankは19であり、非零minorは13759414272*a^12/390625。
P自身のscalar列を追加してもrank19であるから、a≠0の射影PGL orbit次元は18である。
さらにdP/daを追加するとrank20になる。
従ってこの一parameter族のGLとscalarによるsaturationの実像次元はcone20、射影19である。
この次元が他の既知stratumと一致しても、それらへのmembershipは推論しない。

第二sourceのfixed-fiber20方向を実行列g^-1で物理pushbackし、第一sourceと同じPへ揃えた。
両source接空間の和は

    span(C, all9_flag_Lie(P), g^-1_*C)

である。共通Cを消去した150×29行列のQQ(a) rankは17、選択minorは4076863488*a^25/3125。
従って全a≠0で

    dim(T1+T2)=37,
    dim(T1∩T2)=29+29−37=21.

射影ではそれぞれ36と20である。
**37は二つのbranch接空間の和の次元であり、whole imageのembedding dimensionの等号ではない。**
非線形の二branch idealの交わりでは、共通一次式が同じ一次項を持つ実際の共通ideal元へ持ち上がるとは限らない。
また、接空間の一致・次元だけではscheme intersectionのclean性を証明しない。

## 全proper fiberを使った二つのsmooth branch

全CA sourceは、proper flag base上のrank20 vector bundleからの実proper mapである
（`AMBIENT_CA_GLOBAL_FLAG_MODEL_20261103_v63_REVIEW`）。
B18基点P0の**全flag fiber**は一点支持、scheme length2である
（`AMBIENT_B18_QUADRATIC_IMAGE_REDUCTION_20261103_v63_REVIEW`）。
propernessと有限fiberから、P0を含むtarget openで全source mapはfiniteとなる。
finite moduleのfiber lengthの上半連続性により、さらに縮めて全fiber length≤2とできる。

実多項式arcのinverse imageはa=0を含むA1のopenである。
そのa≠0では、既に構成した二つの異なる旗があり、全fiber length≤2なので、
全fiberは**正確に二つの被約点**となる。
両sourceはsmooth29で微分rank29であるから、各completed sourceへのambient形式環の写像はsurjectiveである。
従って二つのkernel I1,I2はsmooth29なformal branchを定める。
cross-tangentの和が37>29なのでI1≠I2である。

ここでglobal birationalityを先に仮定しない。
全finite source代数B1×B2への実completed imageを取ると、reduced CAのcompleted ringは

    A = T/(I1∩I2),  T=ambient formal power series ring,
    T/I1≅C[[29]],  T/I2≅C[[29]].

全fiberを使ったため、第三のformal branchをこの点で落としていない。
射影CAでは二つのsmooth28 branchとなる。
この結論はarc上のa=0を除くあるopen、およびそのG-saturationで成立する。
原点P0のfull completionをこの式に置き換えない。

## 追加の局所帰結と未決部分

二branchのscheme intersection C=T/(I1+I2)のembedding dimensionは21である。
GL×A1の実paired-sourceはrank20で両branchに入り、intersectionにはcone20次元の実族が含まれる。
従ってこの点でのintersectionの次元は20以上21以下だが、全schemeやclean性はまだ決めていない。

この二branch記述が成立する点では

    0→A→C[[29]]⊕C[[29]]→C→0

とdepth lemmaによりdepth A=depth C+1≤22<29。
従ってreduced CAは、この射影19次元のG-saturated self-pair族の一般点でnonnormal・singular・non-CMである。
射影のdepth上限は21である。
これは原点のfull odd idealや全bad locusを分類する主張ではない。
normalization/conductor全体、intersectionの正確な型、R19やNへのmembershipは未決である。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/ambient_CA_origin_self_pair_geometry_20261103_v63.py --output <new-output-directory>

全QQ[a]係数行列、rank witnessとminor:
`reports/AMBIENT_CA_ORIGIN_SELF_PAIR_GEOMETRY_20261103_v63/`。
実行約1秒、制限120秒・2.2GB、GBなし。
