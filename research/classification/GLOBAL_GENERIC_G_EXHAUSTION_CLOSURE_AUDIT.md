# 一般 β に対する G の全分岐閉鎖：大域仮定と次数表の監査

2026-11-08、担当 `/root/plane_collision_review`、PRIMARY G。
本稿は受理済み幾何証明・別作者 receipt の**接続監査**であり、旧計算の再実行ではない。
新 CAS、minor、GB、jet、source probe は実行していない。
全134式と有限 gauge の既存証明を入力とし、固定点・固定 normal form の結論を
全 projective source の結論へ無断で拡大していないかを点検した。

**監査結果：全分岐の接続を PASS とする。** corrected-minus triple の最後の正式 receipt
`PLANE_NONREDUCED_M3_TRIPLE_CORRECTED_G_SEPARATE_REVIEW.md`
（SHA256 `49ab6835d5f9003b70c63b35937e141d2e8c7f7c3b000de889f36552d526a4ae`）も
全三算術 phase・全 source 包含・二つの弱上界壁を含め読解した。
最後の documentary gate は解消済みであり、G の採録に未受理の数学枝を残さない。
旧 plus-sign 九十源の whole coverage は依存に戻さない。

## 1. 正確な対象と、矛盾させる成分の次元

元の S61 全134式の射影を `S61→W=A16`、一般 β の係数体を `K=C(W)` とする。
既存の exact splitter が分離した既知五点を E5、残余代数を B と書く。
対象は `Spec(B⊗K Kbar)` の空性、すなわち G である。

`DOMINANT_COMPONENT_GEOMETRY.md` §1 と
`C0_S61_SLICE_SCALAR_AUDIT_20261106.md` §1–4 によれば、
残余の幾何点があれば、T01/T39/T40 以外の W 優越 full Poisson component C が存在し、

\[
 \dim C=(16+d)+5+4=25+d,\qquad d\in\{0,1,2\}.
\]

ここは PGL5 商前の**全射影 Poisson 空間**の次元である。
5 は固定点での非退化射影評価値、4 は有限 perspective gauge の方向。
55 消去・20 Bianchi の全恒等式と可逆有限 gauge を使い、形式近傍だけの等式ではない。
CSp4 と perspective 群の連結性が成分を保つため、chosen slice だけの次元ではない。
一般 fiber 上界 d≤2 は `research/notes/S61_AFFINE_SURFACE_BOUND_20261102.md` の既存証明を用いる。
「C が W を優越する」は、この rank 四評価・primitive gauge chart の S61 成分が優越する意味であり、
全 PGL5 軌道上に β が一つの大域函数として定義されるとは仮定しない。

ここで「追加成分」は**T01/T39/T40 以外の優越成分**という意味である。
`EXCEPTIONAL_BASE.md` の歴史的52件の一覧を現在の全126件の監査と読み替えない。
元の MPS40 のうち W を優越する三成分と既知五点の対応には同資料の全一次像監査を使う。
後から見つかった成分については、以下の全源次元又は全 primitive β 像の非優越を使う。
新しい成分も含め「既知だから一般 β には現れない」とは仮定しない。

以下、像の上界≤24は常に、plane・surface・mark・normal frame・tensor amplitude を戻し、
実際の同一像 fiber を引いた**大域射影像の上界**を意味する。
固定 source の線形核や一つの tangent rank をこの上界に代用しない。

## 2. 実退化五次 D の全因子型

S61 の基点で rank 四なので、実 Pfaffian 五次 D は非零である。
以下の分類は Poisson scheme の非被約性ではなく、**実際の D の因子と重複度**について行う。

| D の全分岐 | 全源について使う受理済み結論 | 正確な入口 |
|---|---|---|
| 正規 | 全点 T39 | `papers/normal_actual_divisor_recognition_ja.tex`、`DOMINANT_COMPONENT_GEOMETRY.md` §2 の既存依存監査 |
| 被約可約 1+4 | 既知閉和又は大域像≤24 | `HYPERPLANE_QUARTIC_COFACTOR_ENTRANCE.md`、全 nil 十型の SHORT/LONG/COTANGENT/MINIMAL 等の独立 receipt |
| 被約可約 1+1+3 | 既知閉和又は大域像≤24 | `TWO_HYPERPLANES_CUBIC_GENERIC_NONDOMINANCE.md` |
| 被約可約 1+2+2 | 既知閉和又は大域像≤24 | `TWO_QUADRICS_SUPPORT.md`、NONZERO_PAIR/JOINT_SS/PURE_NILPOTENT 全源 |
| 被約可約 2+3 | MPS40 と六つの23次元成分の閉和 | `DEGENERATE_QUADRIC_Q3_AZERO_SUPPORT.md`、Q3_ANIL3、Q4_NIL、SMOOTH_QUADRIC 全源 |
| 被約可約 1+1+1+2 | 全配置が既知成分 | THREE_HYPERPLANES_QUADRIC の dependent/independent 全源と独立 review |
| 被約可約 1^5 | 全配置が T01 又は rank 四不存在 | `FIVE_HYPERPLANE_SUPPORT_REVIEW.md` |
| 非被約、重複 H あり | MPS 閉和又は primitive β 像が proper | `REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` §§3,9,13–16 |
| 非被約、重複 H なし | H+2Q、Q 既約の全点既知被覆 | 同資料 §1、DOUBLE_QUADRIC/RANK3 全支持 |
| 被約既約・非正規 | §3–6 の全曲面次数表 | 以下 |

次数分割は五の全分割で尽くされる。非被約なら重複 prime の次数 e は 2e≤5 から一又は二。
重複 H が無ければ正確に H+2Q であり、複素数体上の既約 Q は rank 三以上。
2+3 の rank 三二十 cofactor 型も有限還元だけで止めていない。
受理済み `R4_QUADRIC_CUBIC_FULL_SUPPORT` は、A0 七型を六23成分へ、Anil3 五型を MPS へ、
残る全型を MPS 又は空へ送り、rank 五・四の全支持も含む。

重複 H の β 議論では、H は rank 四基点を避ける。
その基点を固定し接写像が恒等の perspective で H を無限遠に移す操作は、primitive β を保つ。
全 `J+L+Q`、`Pf(Q)=0` incidence と元 S61 の **β 像の等式**が同資料 §3 で証明されている。
三次尾部を零とした小切片に限定した議論ではない。
全 modular A は可逆、rank 三、rank 二、rank 一、A=0 に分かれ、各零固有値・定数部の全 affine 型を保持する。
rank 二 nil 全五型、rank 一全型、A=0/c非零、A=c=0 の全六葉層閉包まで処理済みである。
最後の E(3)・対数極限は `UNIMODULAR_EXCEPTIONAL_NILPOTENT_INCIDENCE.md` と
`UNIMODULAR_LOGARITHMIC_POLE_BOUNDARY_REVIEW.md` を用い、generic な最高部表示を全境界と同一視しない。
MPS 被覆枝まで含めた全 incidence 像≤15という**強い未証明命題は不要**であり、ここでは使わない。

未知 C は MPS の SNC-meeting component 定理と相異なる成分の極大性により、
MPS40 の有限閉和には含まれない。正規 D の枝も同じ理由で除外される。
従って残りは被約既約・非正規 D だけである。

## 3. 非正規 D が供給する曲面：次数一から六の外は無い

`IRREDUCIBLE_NONNORMAL_QUINTIC_ENTRANCE.md` §§1–2 の全入口を使う。
正規化 ν:X→D、L=ν*O(1)、M=ν*OX/OD の二次元支持成分を S_i、generic length を δ_i>0 とすると、
一般平面切断の正規化種数 g に対し

\[
 \sum_i\delta_i\deg S_i=6-g,\quad 0\leq g\leq5,
 \qquad K_X+C_X\sim0,\quad L^3=5.
\]

D は hypersurface なので S2。非正規性は R1 の破れを含むので、少なくとも一つの integral surface S_i がある。
下位次元の支持を選んでこの表を逃れることは出来ない。
Hamiltonian 導分は正規化・導手を保ち、各 reduced integral S_i は Poisson invariant。
その一般点の restriction rank は二又は零である。
一つでも S_i を選べばよく、全曲面成分が同時に同じ rank であるとは仮定しない。
rank 二では局所 tangent 二次 block が可逆、mixed block は I_S に入る。
実 Pfaffian が I_S² に入るので、その normal 二次 block の一次項も零となり、
intrinsic conormal bracket は可換である。
rank 零では同じ結論に全重複度の modular 補題を使い、rank 二の可逆 block 議論を流用しない。

各 e=deg S_i≤6。e五なら δ_i=1、g零又は一。
g零では残余は degree 一・length 一の plane、g一なら他の二次元導手支持は無い。
e六なら δ_i=1、g零で、他の二次元支持は無い。
この事実が degree 五・六の normalization 分類への**必要仮定**を供給する。
e≤4について g零/一を仮定しない。

S が hyperplane に入り e≥3 なら `I_S=(ell,q_e)` は prime CI。
`I_S^(2)=I_S^2`、D∈I_S^(2) と degree 五から
`F=ell² A3+ell q_e B_(4−e)` となり ell が F を割る。
従ってこの枝は全 restriction rank・全 singular S で空。
退化した span の残りは plane と rank 三／四の integral quadric だけである。

また非plane S の generic transverse multiplicity は正確に二。
三以上なら一般二点の secant line は五次 D に入り、Sec(S) は span P3 又は P4 全体となって矛盾する。
`RANK_ZERO_SURFACE_GENERIC_FIRST_JET.md` §1 の join 証明は rank 非依存である。
**二重度二は node を意味しない**。平方 normal cone と cusp を以後の全表に残す。

## 4. 全非plane S の degree/rank 表

表中の上界は全 PGL5 飽和像。normal、smooth、ACM、clean、ordinary node、vertical 無し、
normalization-flat といった追加仮定を表の全域へ持ち込まない。

| e と全支持 | rank 零 | rank 二 | 全域の独立依存 |
|---|---|---|---|
| 2、smooth quadric | actual H 因子、空 | 同じ全 H 因子、全 R で空 | `SMOOTH_QUADRIC_CONORMAL_FACTOR_SEPARATE_REVIEW.md` |
| 2、quadric cone | 全像≤23 | 全像≤24 | `QUADRIC_CONE_ZERO_REVIEW.md`、`QUADRIC_CONE_RANK_TWO_SEPARATE_REVIEW.md` |
| 3、smooth cubic scroll | 全像≤23 | 全像≤24 | CUBIC_SCROLL_ABELIAN_CONORMAL、A25_BOUNDARY、COORDINATE_AND_POLE_ERRATUM と各独立 review |
| 3、twisted cubic cone | 全像≤24 | 全像≤24 | CUBIC_CONE_ABELIAN_CONORMAL、CUBIC_CONE_RANK_TWO の全別作者 review |
| 4、全 integral 非退化曲面 | 全像≤24、非CI枝は空 | 全像≤24、normal CI22等は空 | `DEGREE_FOUR_RANK_ZERO_SURFACE_EXHAUSTION.md`、regular/singular CI22 reviews、QUARTIC_SURFACES_ON_QUADRIC、NORMAL/NONNORMAL_CI22 reviews |
| 5、g零・残余 plane あり | 全像≤24 | 全像≤24 | RANK_ZERO/RANK_TWO_QUINTIC_SURFACE_GENUS_ZERO の別作者 reviews |
| 5、g一・elliptic normalization | T40 又は全像≤24 | 全空 | RANK_ZERO_GENUS_ONE と RANK_TWO_GENUS_ONE_ELLIPTIC の別作者 reviews |
| 5、g一・canonical index 二 normalization | 全空 | node全像≤24、cusp全空 | rank零 GENUS_ONE、CANONICAL_INDEX_TWO/CUSP_EMPTY の別作者 reviews |
| 6、g零・全 normalization | 全空 | 全像≤24 | `RANK_ZERO_SEXTIC_CONDUCTOR_SEPARATE_REVIEW.md`、`RANK_TWO_SEXTIC_SURFACE_SEPARATE_REVIEW.md` |

e三では minimal-degree 分類により smooth scroll／cone の全二型。
rank零 cone の旧 A1 semisimple 誤ラベルは使わず、corrected V2 の E1 代表を含む三全 strata を使う。
rank二 smooth scroll は canonical roots と残 A25 全二型、四つの full P2 fiber を含む。
同 A25 の incidence 座標・可約 Ω の pole 位置には明示 erratum と別作者 receipt を適用する。

e四の almost-minimal-degree 分類は全外点 projection centre を対象とする。
normal・linearly normal は ACM→CI22、他は Veronese/S13/S22/S04 全射影。
CI22 の pencil は regular 全 Segre 型と all-member singular 全 Kronecker 型を共に含む。
rank二 nonCI scroll/cone の rank 四 containing quadric は全 centre/boundary で確認済みであり、
rank 三 quadric に Weil-class 議論を流用しない。

e五・六では D の正規化 X を分類するのであり、曲面 S の smoothness や normality は不要。
e五の g一には Fukuma の全 polarized model 橋が、elliptic と canonical index 二を尽くす。
zero summands／small contraction／nonterminal／mark 従属・独立の各補集合はそれぞれの receipt に含まれる。
e六の全五 scroll 型113/122/014/023/005も両 rank で処理済み。
rank零 high-cone 排除を rank二へ流用しない。rank二023の proper conductor locus の次元計算には
全 framed projections と actual Aut14 を戻した theorem を使う。
旧 normal-sextic 部分還元に残っていた irregular ruled 枝や未読 Halanay は、この全正規化証明の未処理枝ではない。

## 5. e一・rank零 plane：全 normal order を一括して戻す

`PLANE_RANK_ZERO_WHOLE_SUPPORT_NONDOMINANCE.md` と WHOLE_SUPPORT／DEPENDENT_ORE の
別作者 receipt により、全源は MPS40 又は大域像≤24へ入る。
被約 D の modular flow による conormal 可換性は **全 m** に有効。
従って旧 m二35/30/9 source を使えるが、m二で落とした parameter は必ず復活させる。

traceless co-Higgs 非零の saturation は O(k)、k零又は一。
k一の coefficient rank 三／二は tau零、M(p)零、全 κ 壁を含む high-m 別作者証明で≤24。
rank一と k零は global V 源へ入り、normal ss/nil と [T,V] の全枝を保持する。
leading 非零の m二は既知 SNC 認識又は全壁 bound、高 m は nil b非零、dependent T=tV、
independent/common-image の全 complements により≤24又は rank二以下。
pointwise parallel と constant linear dependence を混同しない。

traceless 部分零では T非零の全十五源が≤23、T零/K非零は全 rank二以下。
T=K=0 の **全三十九 higher 源**は Lie algebroid/CE の別作者証明で≤24。
anchor零の H 因子、全 extension／global cocycles を含む。
これで first symbol 零を非零 co-Higgs source の極限だと仮定する必要は無い。

## 6. e一・rank二 plane：全 R、全 m、全 root 分岐

π|S は O_P2(3) の非零 section f なので、零 cubic R は必ず存在する。
全 R は次の九型で尽くされる。
integral 三型 smooth/node/cusp、reduced reducible 四型 transverse conic+line／tangent conic+line／
triangle／concurrent 三線、nonreduced 二型 2L+M／3L。
node/cusp の非中心係数体と、root-cover の非分裂を初めから split と仮定しない。

plane 沿い exact order は m≥2。m五では F は normal 二変数だけの binary quintic となって分解し、
m六以上は F=0。従って既約被約 D には **m二・三・四だけ**。
`PLANE_IRREDUCIBLE_CUBIC_RANK_TWO_REVIEW.md` は全 integral R を全 normal cone 型で空とする。
残りの六 R に対する表は次である。

| m／normal leading の generic partition | reduced reducible 四 R | nonreduced 二 R |
|---|---|---|
| 2、SF=11 | 全 MPS 又は空 | R2L+M 全空、R3L 全既知 essential三又は空 |
| 2、double=2 | 全像≤23 | moving-root全空、constant-root全像≤23 |
| 3、SF=111 | 非concurrentは全像≤23、concurrentは全像≤24 | 全像≤24又は essential三／空 |
| 3、21 | 同じ全 nonconcurrent／concurrent theorem | 全 generic21 像≤24 |
| 3、triple=3 | 同じ全 nonconcurrent／concurrent theorem | corrected-minus 全 triple 像≤24、二壁は≤23 |
| 4、SF／211／22／31／4 | 全像≤24 | essential三八23成分又は空 |
| 5以上 | 全空 | 全空 |

**m二。** SF は Δ(G2)=c f²、UFD primitive roots から全六 flags 00/01/02/03/11/12。
root basepoints、unused R 成分、全 left-U matrix parameter を残す。
reduced R の whole SF coverage、nonred 全 pole budget、R2L+M SF別作者全空、
R3L SF00/01/02 全 receipt がある。
double は primitive degree零／一、全 irregular scalar trace を残し、
reduced 全四型、nonred moving-root、R3L essential三、最後の R2L+M constant-double 全34源を尽くす。
最終独立 `PLANE_DOUBLE_LINE_M2_CONSTANT_DOUBLE_SEPARATE_REVIEW.md` は
全51壁と critical p=−2/7 の元 row34=−7m²/9まで戻す。
これが受理済み `G_PLANE_RANK_TWO_M2_ALL_NONDOMINANCE` の全域である。

**m三・reduced R。** `PLANE_REDUCED_NONCONCURRENT_M3_ROOT_COVERS_SEPARATE_REVIEW.md` と
NINE_SEPARATE_REVIEW は、transverse 全空、tangent 四・triangle 五の全 higher を処理する。
concurrent は `PLANE_CONCURRENT_M3_WHOLE_NONDOMINANCE.md` と NONORE/NONZERO_QUARTIC/
ZERO_QUARTIC の全別作者 receipt により、二 actual PGL5 orbit と既 essential三全源へ入る。
SF の split/C2/C3/S3、c-dependent／c-independent、21全 flags、triple全 free-r、
K非零211/22/31/4 と K零 SS/Jordan/constant/proportional/central を保持する。
330候補の実列挙を実行済みとはせず、全 strict Q 経路により不要とした証明を使う。
二 actual global orbit に plane六や mark を追加しない。

**m三・nonred R。** SF の全十九 first symbols と全四 Ore modules、c-dependent g二の弱上界24を使う。
21 は `2d_double+d_simple≤2` により 00/01/02/10 の全四 flags。
CONSTANT／MOVING_DOUBLE／MOVING_SIMPLE の別作者 receipt が全非log／全 parameter 壁を含め≤24を与える。
triple は primitive degree零（3d≤2）。ここだけ自由 off-diagonal r の正しい curvature
`dr−(q−p)∧r=0` を使う corrected V2 全九十源を必要とする。
旧 plus-sign whole cover を literal の算術成功で正当化してはいけない。
corrected 全接続→元150/350・30/9・110壁→valid literal全逆橋→全十八 modules／34壁→
唯一新fixedの元 Jacobi ideal の exact z因子、という独立三 phase を用いる。
`PLANE_M3_TRIPLE_VALID_ORE_WEAK_WALL_ERRATUM.md` と addendum は常に併用する。
E1DX_1 の τ=4/21,1/12 は consistent であり、空／proper factor は未証明。
全 source 上界 B≤8、全 Q四、W核六、plane六−projective一から**全像≤23**で充分である。
安全な B≤9でも≤24。E未知を空性に強めない。

**m四。** `PLANE_RANK_TWO_MULTIPLICITY_FOUR_NONDOMINANCE.md` は全 R と
SF／211／22／31／4、image-rank一・高 rank、I/J special、collision、whole Jordan tau零壁まで尽くす。
nonred は NONREDUCED_MULTIPLICITY_FOUR_EIGHT_COMPONENT_COVERAGE、
reduced は scalar/special/repeated の各別作者 theorem を使う。
generic j 又は squarefree のみを全 m四と読み替えていない。

## 7. 有限閉和から G への最後の矢印

各表の全源上界は、必要線形空間の低次 source 数に surface moduli を足し忘れたものではない。
全 source／actual normalizer inverse／その全 rank-jump fiber を保持した後の大域像上界である。
全 parameter family と全複素壁の有限被覆を使い、generic kernel の特殊化だけで済ませない。
有限型像の Zariski 閉包を取っても dimension≤24は変わらない。
proper β 像の枝では、各点・各固定固有値だけの非優越ではなく、全 parameter 源の閉像が proper という受理済み theorem を使う。
この場合は §1 の primitive S61 chart で優越性を比較する。
全 essential三 source の枝には、全八成分認識より弱い
`ESSENTIAL_THREE_QUINTIC_SAITO_BOUND.md` とその別作者 review の大域像≤23だけでも充分である。
八成分の個別相違・新成分性を再計算する義務は G の最後の矢印に含まれない。

したがって C の一般点は、有限個の次の閉集合のいずれかに入らねばならない。
既知 component の閉和、projective dimension≤24の閉像、又は W 像 proper の閉部分。
第一の場合、C の極大既約成分性と既知五点／MPS の一次像に反する。
第二の場合、dim C≥25に反する。
第三の場合、C の W 優越性に反する。
surface を一点に固定して証明してから Hilbert scheme の未知次元を捨てたのではなく、
各表の theorem が全 surface／frame を含む像を抑えていることが本監査の要点である。

よって、最後の corrected triple 全独立 receipt の受理を含め、
**T01/T39/T40 以外の W 優越 component は無く、一般 β の未知幾何支持は空**となる。
既存の exact product `A_K=E5×B` を併用すれば、B が非零の有限型 K 代数なら Kbar 上の幾何点を持つため、
この空性から B=0 も従う。これは既知五点の間に未知の厚みを新しく仮定しない exact splitter の帰結である。
具体的な全底分母又は delta_E の新しい展開・ideal membership 証明書を構成したとはしない。

## 8. 結論の境界と保存

G は一般 β の主張であり、特殊 β 上の全点 E 被覆、全 reduced component の総数、
全 orbit 分類又は全研究完了を意味しない。
consistent な弱上界壁、既知境界に入るか不明な low-dimensional source、
新 componenthood 未証明の actual orbit はそのまま E に残る。
既受理126という下界を上界へ変更しない。

CURRENT/PLAN に歴史的な未完了文が残る場合でも、本稿はその文を削除して証明したことにせず、
後続の exact receipt／明示 erratum／元 hypothesis を優先して接続する。
初回失敗、撤回された whole cover、coarse false report は不変のまま。
依存の exact SHA、最後の receipt の frozen 状態、job gate の scope は同名 report／manifest に保存する。
gate は入力契約の検査であり、数学の独立証明には数えない。
