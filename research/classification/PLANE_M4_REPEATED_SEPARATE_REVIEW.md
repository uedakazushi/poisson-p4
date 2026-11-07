# rank 二 plane・m4 repeated-leading 全支持の独立監査

2026-11-07。主義務 G。別作者による全幾何の読解と、独立 Fraction 実装による
全係数再構成は **PASS**。producer/CAS は import も実行もしていない。
主原稿・global 台帳・Git・旧証明書は変更しない。

## 1. 結論と正確な補集合

rank 四 Poisson 構造の実退化因子 D を既約被約 quintic、S を restriction rank 二の
plane とする。S 沿い位数が四、leading binary quartic が generic repeated root を
持つ全源について、次を独立受理する。

- R が reduced reducible cubic なら、全 PGL5 飽和射影像は二十三次元以下。
  image-rank 一と generic partitions 22/4 は既済の essential-three 八成分へ入る。
  generic211 の高-rank 全六源は像21/22/22以下又は rank 四空、generic31 は下記の
  全 actual D 二正常形で空又は像21以下である。
- generic31 の image-rank 二に対する空性/像21上界は **R が非被約でも成立**する。
  connection の log 性・一意性、R の種類を仮定しない。
- 全 m4 image-rank 一、generic22/4 の essential-three 被覆も R の種類を仮定しない。

既済の squarefree m4 全枝と integral reduced R の空性を合わせれば、
**rank 二 plane・m4・reduced R は全て G 非優越**となる。
本結果の後も、rank 二 plane の m2/m3、及び **nonreduced R の generic211 高-rank**
は未完了である。これらを rank 零 plane の既済 m2 theorem と混同しない。
他支持・全 G・全 E・新成分性は本監査の結論ではない。

有限契約 `jobs/plane_m4_repeated_review.json` は実行前 gate PASS。
独立 code `code/classification_plane_m4_repeated_review.py` は、既済の独立多項式
ライブラリ `classification_plane_m4_special_review.py` だけを再使用する。
結果は `reports/CLASSIFICATION_PLANE_M4_REPEATED_REVIEW_20261107.json`。
全八源の再構成は 6.376 秒、最大 RSS 44,096 KiB である。

## 2. 全 binary quartic 線形空間の fixed multiple root

実 D は normal u,v に関して

```
D = G4(x,y,z;u,v)+g5(u,v),
G4 in H0(OP2(1)) tensor Sym4(C2)
```

である。coefficient image W が一次元なら `G4=L g4` なので D は L,u,v の三変数
だけに依存する。既約 homogeneous quintic は複素二変数以下ではあり得ず、D は
exactly essential-three。既済の全 source theorem により八つの閉23成分に入る。
特定磁気係数、R、partition を仮定しない。

dim W>=2 とし h=gcd(W)、W=hW' と書く。W' は basepoint free であり、一般の二元
a,b を共通零点無しに選べる。標数零の非定数有理関数 a/b の branch values は有限。
一般の a+tb は reduced zeros を持ち、h の有限 root を避ける。
従って W の一般元の multiple roots は h の固定 multiple roots に限る。
これにより generic22/4 なら deg h=4、dim W=1。generic31 なら fixed triple root、
generic211 なら fixed double root を持つ。これは P1 上の elementary Bertini であり、
singular quartic pencil の未読分類に依存しない。

## 3. generic211 の全六正常形と接続一意性

fixed double root を u として
`G4=u²(Au²+2Buv+Cv²)`、`T=AC-B²` と書く。
generic211 の条件は C,T 非零である。U=S\R 上で G4/f は flat conormal connection
の parallel tensor。Hamiltonian flows が Pi² と PiS を保存する associated-graded
恒等式に基づく既済入口で、simple-pole/log/ordinary-node を仮定しない。
非零 parallel tensor の root partition は U 上で変わらないので、C=0 と T=0 の
全 collision divisor は R に支持される。basepoint で parallel tensor が零なら U
全体で零になり、generic leading 非零に反する。unused R components は許す。

rank(A,B,C)=3 なら T は irreducible nondegenerate conic、C=0 はその tangent。
reduced cubic R には両方が必要なので tangent conic+line だけである。
rank 二では regular symmetric2 pencil の congruence normal forms は、distinct
determinant roots の `Lr²+Ms²` と、Jordan の `Lr²+2Mrs` の二型に尽くされる。
invertible symmetric 元を選び self-adjoint endomorphism を Jordan 化する証明で、
Jordan 境界を捨てない。determinant 恒等零では一般 Q が square となり211ではない。

第一型の marked root l が r 又は s なら used lines は二つ、どちらとも異なれば
正規化して l=r+s、used lines は L,M,L+M の三 concurrent lines。
第二型で l=r は31へ送り、l≠r は upper-triangular stabilizer により l=s。
used lines は L,M の二つである。従って reduced reducible R の全高-rank 源は
次の六型に尽くされる。

| R | Q の全正常形（fixed double root は u） | 全 parallel kernel 次元 |
|---|---|---:|
| xyz | `x u²+y v²` | 2 |
| xyz | `y uv+x v²` | 2 |
| xy(x-y) | `x u²+y v²` | 2 |
| xy(x-y) | `y uv+x v²` | 2 |
| xy(x-y) | `x u²-2xuv+(x-y)v²` | 1 |
| z(xz-y²) | `x u²+2yuv+zv²` | 1 |

全 leading は u²Q。最初の四型の二振幅、後二型の共通振幅は全 kernel に保持する。
振幅零は generic partition で別枝に入り、勝手に捨てない。finite line permutations、
normal GL2、base changes は actual invertible maps であり全 chart を覆う。
transverse conic+line にはこの高-rank211源は無い。

三つの distinct roots を持つ binary quartic に対する twisted GL2 infinitesimal
stabilizer は零。projective 変換は三 root を固定して零、残る scalar は
Sym4 tensor det^-1 上で非零 weight 二を持つためである。従って leading/f から
全 rational connection が一意に復元され、nonlog/off-diagonal 枝も含む。
normal covectors 上の convention を固定すると

```
beta_uu = (1/2)dlog f-(1/4)dlog T,
beta_uv = 0,
beta_vv = (1/2)dlog f-dlog C+(1/4)dlog T,
beta_vu = ((B/2)dlog T-dB)/C.
```

全六型で f beta は polynomial one-form となる。Jordan 型で f/T 自体が
polynomial であるとは仮定しない。独立実装では上式を f,C,T で cross-multiply
した恒等式として全係数検証し、全十二 curvature numerator `d beta-beta wedge beta`
の零も確認した。旧 producer の f/T 実装失敗は成功とは扱わない。

## 4. 全 P0 と全30/9源の再構成

`tau=beta-(dlog f/3)I` とし、tracefree linear base fields から
`iE iV vol=3f tau` を解く。完全18×8写像は rank 八であり、保存された全四 entry の
逆解を再検証した。base Jacobian Pf と、全 raw mixed brackets
`{xi,na}=-sum_b nb(Vab)i` を組み直した。
raw normal-normal quadratic K の全三係数は実 D2=0 により一意に決まり、
全126係数の rank 三 forcing matrix/minor と全代入を確認した。
triangle diagonal は K=uv/4、triangle Jordan は K=-uv/2、tangent net は K=3u²/4。
concurrent 三型だけが K=0 である。入口で K を零にしていない。
全150 P0、全350 Jacobi 零、実 D 零を再構成した。

leading parallel equation の全十五係数、全 matrix/kernel/minor と annihilator も
再構成し、前表の全振幅を証明した。higher tensor の weight 一は raw36/Euler6 の
完全 quotient30、weight 二は完全9であることを独立 source 関数で再確認した。
全 original350 weight一 Jacobi と D3 の全24係数を合わせた374×30 matrix、
weight二の全350と全 leading annihilator shape rows、全 kernel/minor を再計算した。
さらに全 P1 quadratic polarization と元 row の係数を照合した。

P1 が全 K1 に入ると、weight 二の式は同じ定数 operator L による
`L(P2)+Q(P1)=0`。各 fiber は空又は ker L の affine translate であり、
generic P1 の rank を特殊 P1 へ外挿していない。全 source 上一様に
dim source<=dim K1+dim ker L である。

| 全型 | dim K1 | dim ker L | fixed P0 source上界 | actual 群次元 | 全像上界 |
|---|---:|---:|---:|---:|---:|
| triangle diagonal | 6 | 1 | 7 | 10 | 21 |
| triangle Jordan | 6 | 2 | 8 | 10 | 22 |
| tangent conic net | 6 | 1 | 7 | 9 | 22 |
| concurrent diagonal/Jordan | 12 | 6 | 18 | 11 | rank四空 |
| concurrent three-used | 12 | 6 | 18 | 10 | rank四空 |

concurrent 三型では P0、全十二 K1、全六 ker L は全て
`partial_z wedge W`、`partial_z W=0` である。Q は全元 row で恒等零なので
任意の P1 に対して P2 は全 ker L に入り、全十八 higher 係数を戻した全源が
rank<=2 となる。粗上界31/31/32を成功と呼ばず、この全構造で rank 四空性を得る。

## 5. 実群の fiber と振幅境界

triangle の block-diagonal 実群は base torus二＋normal diagonal二。
normal-to-base shear六を加え十次元であり、connection を projective に保つ。
tangent net では marked-point Borel の二次元作用を base symmetric-square と
normal 作用へ**同時**に持ち上げ、normal scalar一と shear六を加え九次元。
base と normal の誘導作用を独立に数えない。concurrent の全三型についても
保存された実群の11/11/10はこの同時作用によるが、空性にはその減算は不要である。

全 group は fixed-P0 全 Poisson source X を保つ。block diagonal は leading と
connection を、shear は associated graded P0 を保ち、変化は完全 P1/P2 内にある。
raw K は実 D2=0 の一意解なので群変換後も同じ P0 へ戻る。全 leading kernel の
振幅境界も保持される。PGL5×X への右作用は PGL5 frame 因子から group 元が
復元されるので自由であり、像次元は24+dim X-dim group以下。
個別 Poisson tensor の stabilizer が零である必要はない。
これで上表の全像上界が全支持に適用される。

## 6. generic31 全 actual D の二正常形

本節は nonreduced R を含む。fixed triple root と leading image-rank 二から

```
D=u³(Au+Bv)+g5(u,v),  A,B independent.
```

κ=g5(0,1) は非零。零なら D は u で割れ、既約性に反する。
A,B に normal linear forms を足す actual invertible shear で u-degree>=3 の
三項を除き、`alpha u²v³+beta uv⁴+kappa v⁵` を残す。
v の shear `v -> v-beta u/(5kappa)` と同じ A,B shear の再適用で uv⁴ を消す。
新しい u²v³ 係数は `alpha-2beta²/(5kappa)`。
非零 scalar、A,B の共通尺度で κ=1、u の非零尺度で残る係数を ε=0又は1へ直す。
必要な平方根は C 上の有限選択であり、零/非零両枝を保持する。
従って全 actual D は

```
F_eps=u⁴A+u³vB+eps u²v³+v⁵,  eps=0,1
```

へ尽くされる。残る第五座標を w とする。両 F は A に関する primitive linear
polynomial なので絶対既約、partial derivatives の独立性により essential-four。
ここで二正常形の相互相違や Poisson 構造の存在は主張しない。

## 7. F0 は全 source で rank 二以下

quadratic gradient syzygy を δ=P∂A+Q∂B+R∂u+S∂v とする。
F0 の gradient 式で u=0 とすると S=uL、再度 u=0 とすることで S=c u²。
さらに割って u=0 とすれば `3B R|u=0+5c v³=0`。B と v は互いに素なので c=0、
R=uM。残る式は

```
u(P+4AM)+v(Q+3BM)=0.
```

従って `P+4AM=vN, Q+3BM=-uN`、M,N は linear であり、全 AR2 は

```
S1 Z + S1 Y,
Z=v partial_A-u partial_B,
Y=-4A partial_A-3B partial_B+u partial_u
```

に等しい。AR0=0、AR1=span(Z,Y) である。独立実装は全 degree0/1/2 gradient
matrix、完全 kernel、非零 rank minors を再計算し、次元0/2/8と八つの表示方向の
独立性・全代入零を確認した。点での image rank の推測ではない。

canonical divergence-zero strict lift は既済の Pfaffian-Casimir identity
`Pi#dF=0` を満たす。この identity と lift の GL5 equivariance は essential-variable
数に依存せず、本 essential-four F に適用できる。
各 row の base 四成分を w の次数で展開すると AR2、AR1、AR0 に入る。
従って Pi の image を base へ射影したものは pointwise span(Z,Y) 内。
∂w を戻しても image 次元<=3、skew rank は偶数なので全 rank<=2 となる。
full nonlinear Jacobi、R の被約性、connection 一意性を仮定しない。

## 8. F1 は全150 strict Casimir 源から像二十一以下

各 ε で全150 quadratic skew 係数を保ち、`Pi#dF_eps=0` の degree六全係数
5×210と、div Pi=0 の全 linear 係数5×5を合わせた**全1075×150** matrixを
独立に再構成した。全 row/column labels、全 matrix、rank minor、完全 kernel の
独立性と全積零を検証した。F0 の kernel は五、F1 は四次元である。
canonical lift を取った後なので、この150源から Euler24を再度引くことはしない。
Jacobi・actual Pfaffian 条件は必要源を狭めるだけである。

Z は F1 も保存する genuine 非零 unipotent PGL4 作用で、全多項式恒等式
`Z(F_eps)=0` を確認した。従って actual-D PGL4 orbit の次元<=14。
essential-four variables の choice Gr(4,5) の四次元を戻して actual-D family<=18。
canonical lift の equivariance により固定 [D] 上の必要源次元は全 orbit で四。
projectivization の一をここで引き、全 Poisson 像は

```
18+(4-1)=21
```

以下となる。固定点の ambient tangent、chosen slice、generic finite prime の
結果を全支持へ外挿していない。以上で31の全 high-rank leading 源が閉じる。

## 9. 凍結入力と再現

211入口は `inputs/plane_m4_repeated_211_manifest_20261107.json`、
31入口は `inputs/plane_m4_triple_root_quintics_manifest_20261107.json`。
全 manifest の入力 hash を検査した。211 の実行 source SHA256
`ca361245a1915735e31c5bbb8acd1170a0492b0ed73d6ca42f5fb58ae6c6b2fb`
は、後続31modeが append された current runner とは異なる。
旧211 report の immutable `producer_source_bytes` を hash で復元して bind し、
current bytes を旧実行と同一とは扱わない。両 report の全 executed snapshot hash を検証した。

```
python scripts/classification.py gate research/classification/jobs/plane_m4_repeated_review.json
python code/classification_plane_m4_repeated_review.py
```

全 source 幾何、全 original equations、actual normalizer と各補集合を合わせて
独立受理する。gate だけが数学の証明であるとはしない。全 G/E の完了、新成分数、
nonreduced211 や m2/m3 の成功はこの receipt から導かない。
