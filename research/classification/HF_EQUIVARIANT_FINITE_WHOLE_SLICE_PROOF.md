# 二つのHF anchor：normal character−2が全formal sliceを有限graphへ強制する

2026-11-08、PRIMARY E。producer純証明。後続graph算術はまだ実行していない。
K4aの既受理 `K4A_EQUIVARIANT_FINITE_SLICE_PROOF.md` と同じétale/IFT機構を、
今回の二つの固定HF anchorと**完全に有限な全次数**へ適用する。
これはformal jetの長さを増やす実験ではない。

## 1. 全ambientで確定した入力

座標(a,b,c,u,v)。Aはk−1/2,p2/5/ε0のp1=q1=1、Bはk−2,p−1/2/ε0のp5=p2=1、他higher係数零。
元tensorの整数20倍／2倍だけを用いた。
両方とも元350Jacobi零、fullD350×150 rank101、Euler25 rank24、実Euler/GL/scalar M51 rank48。
全行列、全right kernel、free identity、非零最大minor、全DM零、全五Pfaffianを保存した。
射影actual orbitは23次元、射影接空間は24次元。この時点でcomponenthood又は親を主張しない。

両actual点をrawで固定するnormalized torusは

```
A: s^(-1) diag(s²,s^(-1),s³,1,s)_*,
B: s^(-1) diag(s²,s²,s³,1,s)_*.
```

入力coefficient x_ij,e のcharacterは `wi+wj−e·w−1`。
Aの全150character範囲は[−8,6]、Bは[−6,4]。
全350equationのcharacterは `sum(w_quad)−e·w−2`。
同じD/kernel/Mのeigenspacesから計算した唯一のnormal quotientは、両方とも**character−2、一dim**である。
新しいnormal operator又はprofileを追加して得た値ではない。

## 2. 真の102-dimensional sliceを既存certificateだけから作る

Dの既存非零101minorのcoordinate columnsをU、equation rowsをRとする。
Uは101coordinate directionsでありtorus安定、D_R,Uは可逆である。
kerDは49次元、Mのimage Gは48次元。
character−2のkerD vector NをGの外に取る。

実際のQQ certificateには、Mの既存48minorだけでNの非group性を示せる。
その48columns G0と48rows Sを取り、候補N_kに対し

`N_k−G0*(G0_S)^(-1)*(N_k)_S`

を計算する。character−2 kernelに非零residualが必ずあり、その最初の非零coordinateから
全51columnsを殺すlinear functionalとN上の非零値を得る。
新150minorや別rank profileを作らない。これは既存完全kernelのbasis選択とcertificate読戻しである。

D*Uの独立性とkerD=G⊕CNより、全raw150spaceは

`G ⊕ U ⊕ CN`

という直和である。従ってP+U+CNは真の102dim affine sliceで、
旧orbit tangentとの補空間である。Nは抽象normal classだけでなくactual quadratic bivectorとして保持する。

## 3. IFT graphの全次数を最初に固定する

sliceをP+tN+U*f(t)と書く。選択101元式F_Rのimplicit derivative D_R,Uは可逆なので、
completionでは選択式の解は一意graph f_i(t)である。
f_i(0)=0、DN=0からf_i'(0)=0。
sliceも選択式もtorus安定、graphは一意なので、各f_iはtorus-equivariantである。
normal parameterはt→s^(-2)tだから、coefficient characterχ_iのformal functionは

`f_i(t)=constant*t^j, χ_i=−2j`

でしかあり得ない。j零／一は既に零、正character又は奇数negative characterにはformal項がない。
全150characterの最小値から、**全completionのgraphそのもの**が

```
A: P+tN+t²T2+t³T3+t⁴T4,
B: P+tN+t²T2+t³T3
```

に一致する。ここでTjはU内、character−2jの方向だけにある。
それぞれj=2,…,4又は2,…,3を一度だけ、同じD_R,U inverseで解けば全graphが求まる。
計算を見ながらjet次数を延長する余地はない。

元350式はPiの係数についてquadraticなので、このliteral polynomial graph上の全式のt次数は
Aで高々8、Bで高々6である。全350の全0..8／0..6係数を検査すれば、
未確認のhigher termは存在しない。

## 4. 全reduced supportとactual親の二つの結論

全350がliteral graph上で零ならslice completionはC[[t]]であり、実polynomial Poisson曲線が存在する。
非零residualが一つでもあれば、全式をgraphへ制限したidealはC[[t]]内の非零idealなので、
そのreduced supportは原点だけである。必要なのはこの二者択一であり、exactnil指数は求めない。

Euler translations／GL5／nonzero tensor scaleの実群作用を考える。
Mの既存48independent columnsに対応する48oneparameter algebraic group factorsを取り、
そのproductと上のsliceからraw150spaceへのmapを作る。
原点で微分はG⊕U⊕CNの直和により可逆なので、このmapはétale。
Poisson schemeは実群作用で保たれ、completionは48smooth group directionsとsliceのproductである。
closed orbit、reductive stabilizer又はLuna theoremは仮定しない。

全residual零の場合、raw局所schemeはsmooth49、Euler24＋射影scale一を除きsmooth unique24componentを得る。
非零の場合、reduced raw局所supportはsmooth48、射影supportは23で、actual23orbitclosureがfull reduced componentになる。
非被約厚みが残る可能性をsmoothschemeと混同しない。

全residual零のactual curveでは、normalized torusはt≠0の全点を互いにactual PGL同型にする。
s→s²がC*上surjectiveなのでnonzero tは全一軌道である。
étale productのdimension下界と同じcurve-torus作用から、その一般actualorbitは24dimで親componentに稠密となる。
別のgenericpoint行列を追加する必要はない。

## 5. 元の全A10 familyを戻す

`PLANE_M2_HF_SIX_EXCLUSIONS_AND_TWO_WHOLE_FAMILIES.md` の閉image論証により、
各全A10のPGL飽和閉包はそのanchorの23orbitclosureと一致する。
上のどちらの結論でも、このorbitclosureは明示したfull23component又はunique24親の境界へ入り、
全A10 familyの全pointwise E帰属がその同じ閉族へ戻る。
全sourceでのtransitivity、unknownfiberのflatness、genericbetaの仮定は不要である。

後続jobはこの全graph/全350有限判定だけとし、他point/profile/minor/jet系列には進まない。
実24親が出た場合の既MPS component名、又はfull23の場合のprior129相違は統合担当の別義務である。
