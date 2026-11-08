# Eの全点被覆：G完了後の全source inventory

2026-11-08、PRIMARY E、担当 /root/log_manuscript。基準は129成分／1044claimsの受理済みcut。
本稿は既存の全源定理・独立receiptの**結論の射程の一覧**であり、新claim・計算・成分性証明ではない。
主原稿・global台帳・Git・旧証明書を変更せず、数学を再実行していない。

入口は[GLOBAL_GENERIC_G_EXHAUSTION_CLOSURE_AUDIT](GLOBAL_GENERIC_G_EXHAUSTION_CLOSURE_AUDIT.md)、
[93-file manifest](inputs/global_generic_g_exhaustion_closure_audit_manifest_20261108.json)、
[34 accepted-claim projection](inputs/global_generic_g_exhaustion_accepted_branches_20261108.json)。
**93は依存ファイル数であり、source族や未解決課題の数ではない。**
snapshot後のcorrected-triple最終receipt、K22陽性、K4a/T37及びclaim1022–1044の全源認識も併用する。
旧note冒頭の「未了／review前」は後続の受理済みreceiptとerratumを優先して読む。
旧52件のexceptional-base表を現在の全E表とは扱わない。

## 1. 記号と数え方

| 記号 | 全sourceについて証明した結論 | Eに残ること |
|---|---|---|
| K | 名前付き証明済み閉成分の有限和へ全点包含 | 同sourceの帰属を再開しない |
| Ø | 指定する実D・rank・位数の仮定の下で全空、rank≤2又はD可約 | 仮定外を適切な別行へ戻す |
| P | 明示actual軌道閉包が全ambient成分で、全源がそこへ入るがcatalogue名／相違は未了 | 成分の識別・重複排除、必要なら採録 |
| B | 全像≤24、primitive β像proper、必要supersource上界などGに十分な結論だけ | 全実源を既知／新成分へ覆うか、全空性を証明する |

**Bは非空性も意味しない。** 必要線形supersourceしか解いていない場合は、実rank四・既約D点の存在自体も未決定である。
小さいkernel、同一jetの実fiber、余核零、properなβ像はKではない。
全actualDのproper因子という受理済み恒等式は、既約Dの行をØにしてよい。

同じtensorが複数の不変曲面やHを持つため、以下は重複し得る**被覆用全sourceの入口**であり、相異なる成分数ではない。
他のK経路へ全source包含を証明できれば重複義務を消せるが、その橋を仮定しない。
各行は全normalizer逆、反転係数零の補集合、全higher、inhomogeneous particularを含む。

Gは完了し、一般βのexact splitter残余Bは零、一般S61代数は既知E5である。
これは特殊β上の全点membershipや低次元成分を消さない。
**plane三源／stable二源／Ga9 linear-w／nilshortだけがE残余、という一覧は不完全である。**

## 2. 実退化五次Dの最上位全分岐

一般rank≤2の全支持は既六成分・既rank四閉包への既存全支持定理で処理済み。以下は一般rank四。

| D | 全点状況 |
|---|---|
| 正規D | K=T39。[normal actual divisor theorem](../../papers/normal_actual_divisor_recognition_ja.tex)。Dredの正規性へ拡張しない |
| 五相異H | K=T01又はØ。[FIVE_HYPERPLANE_SUPPORT_REVIEW](FIVE_HYPERPLANE_SUPPORT_REVIEW.md) |
| 三相異H＋既約Q | K。全dependent/independent配置、全磁気。[DEPENDENT](THREE_HYPERPLANES_QUADRIC_DEPENDENT.md)／[INDEPENDENT](THREE_HYPERPLANES_QUADRIC_INDEPENDENT.md) |
| 既約Q＋既約C | K=MPS＋既六23成分。[A0](DEGENERATE_QUADRIC_Q3_AZERO_SUPPORT.md)／[Anil3](DEGENERATE_QUADRIC_Q3_ANIL3_SUPPORT.md)と全rank4/5支持 |
| 重複H無しの非被約D | 正確にH＋2Q、Q既約。全rank≥3・全配置K。[DOUBLE_QUADRIC](DOUBLE_QUADRIC_SUPPORT.md)／[RANK_THREE](RANK_THREE_QUADRIC_SUPPORT.md) |
| 重複Hあり | K/B混合。§3の全modular型が必要 |
| 被約H＋既約quartic | 大部分K、残B。§4 |
| 被約二H＋既約cubic | 全非plane K、plane三whole六磁気B。§4 |
| 被約H＋相異二既約Q | plane全十候補K、jointSSrank一・pure-nil factor-specific枝にB。§4 |
| 被約既約非正規quintic | 導手曲面次数1..6・rank零/二のK/Ø/B/P。§5–7 |

## 3. 重複Hの全modular型

全源はh-Casimir liftのP=C+L+Q、Pf(Q)=0、全Jacobi。
symplectic基点／finite perspectiveとの橋は[REMAINING_SUPPORT §3](REMAINING_DEGENERACY_DIVISOR_SUPPORT.md)。
固定評価・Taylor fiberを引くG論法をE帰属としない。

| Z=c+Ax | 状態と全source入口 |
|---|---|
| A可逆 | K。[REMAINING §§5–9](REMAINING_DEGENERACY_DIVISOR_SUPPORT.md)、[MODULAR_SUPPORT](MULTIPLE_HYPERPLANE_MODULAR_SUPPORT.md)。全一次共鳴・対和零・退化C・非半単純を含む |
| rankA3、全非零半単純部 | K。同noteの零点あり/無し、±a＋Jordan2を含む全実MPS変形 |
| rankA3 nilJordan4 | K。[NILPOTENT_SUPPORT](MULTIPLE_HYPERPLANE_NILPOTENT_SUPPORT.md)：全五源→既五22閉包 |
| rankA2半単純 | K。[LOW_RANK §2](MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md)：c全余核・全三共鳴・g定数/非定数/零 |
| nil31、c∈imA | K。[MODULAR_NIL31_COMPONENT](MODULAR_NIL31_COMPONENT.md)：n7開と全閉補集合 |
| nil31、長いtop | K。[MODULAR_NIL31_LONG_SUPPORT](MODULAR_NIL31_LONG_SUPPORT.md)：全六源、最後r=−1も全源変形 |
| nil31、短いtop | **r=2/3全源B**。[SHORT §§7–8](MODULAR_NIL31_SHORT_SUPPORT.md)：元六変数→(t,d)、t≠0、td≠0一actual軌道、d=0を同実源に保持。r=1/3主枝は後続§8でMPS、追加枝T01、他三源は既22。一般点の親／成分性を得れば全源へ閉性を戻せる |
| nil22、c∉imA | K。[MODULAR_NIL22_TOP_SUPPORT](MODULAR_NIL22_TOP_SUPPORT.md)：全八源、全MPS/T19/二22。途中の古い未了文を最終状況としない |
| nil22、c∈imA | K＋**stable二B**。[NIL22_IMAGE_SOURCE_REVIEW](NIL22_IMAGE_SOURCE_REVIEW.md)：全十一Hessian源、GL2全逆、任意F0六係数。half/third等、regular limits、doubleB/C、quadraticA、doubleAの後続陽性を合わせ、残stable_line_conic／stable_three_lines |
| rankA1、c∈imA | K。[LOW_RANK §4](MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md)：b≠0は既nil31成分、b=0はCA/CB |
| rankA1、c∉imA | **K/B混合**。同§5長鎖/短鎖全源。長鎖g定数非零・H3=0はNEW25へK。g=0の自由F二次、非定数g=zの全U,K、定数gの三次共鳴H3等は同一jet fiber／β像上界のみB |
| A0、c≠0 | **B**。同§3／remaining §14。全P=∂t∧Y+Jac(f)、divY=κ、Yf=κf、degY≤2、degf≤3。β(c)=0はmembershipではない |
| A=c=0 | 以下の全symplectic延長でK/Ø/B混合。全点をbeta零／Q45に制限しない |
| 別の横断条件[J,Q]=0 | **B**。remaining §4の全Lie4源＋全Q。固定評価β像≤12だけ。上の型と重複する必要sourceとして全係数を保持 |

### A=c=0：最高部六閉包の全逆像

全問題はclosed二形式α=J+β+γと、J²≠0、J∧β=0、
2J∧γ+β²=0、β∧γ=γ²=0の全式である。最高部六閉包へのG被覆は**低次延長のK定理ではない**。

| 最高部入口 | 現在の全源結論／E義務 |
|---|---|
| γ零／共通因子 | B：全低次β・非退化J。primitive評価の核という必要閉条件だけ |
| R13／linear Casimir | ØとB。[UNIMODULAR_LINEAR_HAMILTONIAN_SUPPORT](UNIMODULAR_LINEAR_HAMILTONIAN_SUPPORT.md)：零点無しaffine Hamiltonian、semi/constantとnil六対、五つの全45係数線形源 |
| R22／quadratic pencil | ØとB。[UNIMODULAR_REGULAR_QUADRATIC_PENCIL_SUPPORT](UNIMODULAR_REGULAR_QUADRATIC_PENCIL_SUPPORT.md)：全Segre/common-kernel/Kronecker、全24一次＋六定数 |
| S23／pullback | B又は前行へ。[REMAINING §16](REMAINING_DEGENERACY_DIVISOR_SUPPORT.md)：w dw等も保持した全低次延長 |
| L1111/L112 | SNC開K、pole境界Ø/B。[LOGARITHMIC_BOUNDARY_REVIEW](UNIMODULAR_LOGARITHMIC_POLE_BOUNDARY_REVIEW.md)：多くのx²yz/x²qは延長無し、混合q=xz+yw全源は評価rank≤3だけ |
| E3閉包 | Ø/Bへ還元。[EXCEPTIONAL_NILPOTENT_INCIDENCE](UNIMODULAR_EXCEPTIONAL_NILPOTENT_INCIDENCE.md)：全四Jordan・denseL112/linear Casimir/pullbackへの全延長。Bへの還元をKとしない |

## 4. 被約可約Dの残る全源

### H＋既約quartic

[COFACTOR_ENTRANCE](HYPERPLANE_QUARTIC_COFACTOR_ENTRANCE.md)の半単純24代表／nil十型と
[全90係数review](HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_REVIEW.md)を使う。
全半単純はMPS又は空。nil1は既22＋T39、nil2/3は空、nil4はMPS、nil5は既22＋T40、
nil7全七源は既成分/MPS、nil8はMINIMAL_M、nil10全translation orbitは半単純実曲線でMPS。
quartic既約性を使った証明をH＋二Qへ流用しない。

nil6 base括弧f非零は[TWO_QUADRICS_NIL22_BASE_SUPPORT](TWO_QUADRICS_NIL22_BASE_SUPPORT.md)の
全二閉源A/Bから重複Hへ戻るか旧H_QUARTIC_ESSENTIAL5へK。
f零は[COTANGENT_G_BOUND](HYPERPLANE_QUARTIC_NIL22_COTANGENT_G_BOUND.md)の全24源、
[MARKED_ENTRANCE](PLANE_HYPERPLANE_QUARTIC_MARKED_ENTRANCE.md)へ入る。
μ6/r2十一とμ5/r1五は[MU6](PLANE_HYPERPLANE_QUARTIC_MU6_SOURCE_SEPARATE_REVIEW.md)／
[MU5](PLANE_HYPERPLANE_QUARTIC_MU5_SOURCE_SEPARATE_REVIEW.md)／
[二normal最終review](MARKED_QUARTIC_TWO_NORMAL_SEPARATE_REVIEW.md)により、全frame・六磁気・zero invariantまでK。
μ6/r1を落とさず、[MU6_MARKED_LIST §§1–4](PLANE_HYPERPLANE_QUARTIC_MU6_MARKED_LIST.md)から
A6 tangent／A4+A2 tangentはGa9power1／power3へ、E6 tangentはτ13空枝又はGa9linear-wへ戻る。
3A2にr1は無い。power1/3は[全六磁気cuspidal24橋](GA9_CUSPIDAL_QUARTIC_COMPONENT.md)からK。
nil9 Ga9の他源も既23/空へ入り、**残Ga9 c0 linear-w全源**がB。
そのdense23軌道のT39/T40等への陰性を陽性被覆とはしない。

### 二H＋既約cubic

[TWO_HYPERPLANES_CUBIC_SUPPORT](TWO_HYPERPLANES_CUBIC_SUPPORT.md)の全非plane十一型はK。
translation二22閉包、j21(h,x)、j21(0,h+x)の全n6補集合、p0→T06の全元源・逆を戻す。
平面十三配置中、未帰属は以下の**全六磁気三source**。

| source名 | E義務 |
|---|---|
| node_flex_tangent_through_flex | whole23閉包の親又は成分性 |
| cusp_flex_secant_through_flex | 同上 |
| cusp_cusp_tangent_cusp_secant | 同上。I≠0全lambdaのglobalT39/T40陰性をI0へ延長しない |

全括弧／Saito frameは[TWO_LINES_PLANE_CUBIC_SUPPORT](TWO_LINES_PLANE_CUBIC_SUPPORT.md)。
旧T36–38 uniform normal-cone陰性は指定入口だけであり、whole親の不存在ではない。

### H＋相異二既約Q

λ非零・jointSSrank二・pair-zero等で実Kとした全sourceと、本質三planeの十候補→七free源はK。
[TWO_QUADRICS_PLANE_WHOLE_COVERAGE](TWO_QUADRICS_PLANE_WHOLE_COVERAGE.md)を次の非plane/別factor源へ拡張しない。

| Bのfactor-specific入口 | fullsource／必要unit |
|---|---|
| λ0 jointSSrank一B1 | [JOINT_SS_RANK1 §3](TWO_QUADRICS_JOINT_SS_RANK1_SUPPORT.md)：六自由n0,n3,n4,n5,n8,n9、n0n5n8≠0 |
| 同B2 | 同§4：n0非零／零の二つの全七parameter源、n6n9≠0 |
| 同B3 | 同§5：n1非零／零の二つの全五parameter源、指定n6/n7等のunit |
| pure-nil型1 regular c0 | [PURE_NIL §2](TWO_QUADRICS_PURE_NILPOTENT_SUPPORT.md)：n12=−1/2,n10=0、n0≠0、n0,n3,n4＋実N² undo n5。Pf=I2²/48−n0w²/2 |
| pure-nil型4 | 同§1／[NIL31_SHORT](HYPERPLANE_QUARTIC_NIL31_SHORT_SUPPORT.md)：全factor-twoQ源≤23。quartic既約のKを全factor枝へ外挿しない |
| pure-nil型7 | [NIL22_LONG_SUPPORT_REVIEW](HYPERPLANE_QUARTIC_NIL22_LONG_SUPPORT_REVIEW.md)：全20potential/670lift/criticalcore、≤24 |
| pure-nil型10 top2 | [PURE_NIL §6](TWO_QUADRICS_PURE_NILPOTENT_SUPPORT.md)：全canonical三parameterδ,a≠0,κ、Pf=δR²/2+az²、≤24 |

A1はØ、A2はplaneKへ戻る。pure-nil2/3/5/8はØ、6/9はplaneK。
残表では相異二既約Qという実factor条件を全supersourceへ戻してから親を認識する。

## 5. 既約非正規quintic：非plane導手曲面の全次数

共通入口は[IRREDUCIBLE_NONNORMAL_QUINTIC_ENTRANCE](IRREDUCIBLE_NONNORMAL_QUINTIC_ENTRANCE.md)。
Σδ_i degS_i=6−gの全支持を保持する。Sのnormality/smoothness、ordinary node、clean、
vertical無しを表の全域に加えない。

| degS／restriction rank | K／Ø部分 | **Bの全source**と次の認識対象 |
|---|---|---|
| hyperplane内deg≥3 | Ø | symbolic-square次数でH因子。補集合無し |
| 2 smoothquadric、両rank | Ø | [CONORMAL_FACTOR_REVIEW](SMOOTH_QUADRIC_CONORMAL_FACTOR_SEPARATE_REVIEW.md)：全R/vertical |
| 2 quadriccone、rank0 | **全Ø** | [ZERO_ACTUAL_FACTOR](QUADRIC_CONE_ZERO_ACTUAL_FACTOR.md)／[root独立読解](QUADRIC_CONE_CI22_AND_TRIPLE_Y_ROOT_REVIEW.md)：全semi/nil・全十五Y・全共鳴のfull affine bad-rowが一次式で割れ、全actualDにH因子。λ零/C零及び全逆も保持。可約D全源の所属とは区別する |
| 2 同rank2 | C0/Fv零等Ø | [RANK_TWO_REVIEW](QUADRIC_CONE_RANK_TWO_SEPARATE_REVIEW.md)：全33源、正しいHamiltonian pencil→D0後だけOre。全κ,s、十二semline＋七nilline、全72+13根・C0全五根。最後nil(0,0)のp1(p0+p2/2)は次元一低下でありKではない |
| 3 smoothscroll、rank0 | **全Ø** | [修正V2全支持と独立review](CUBIC_SCROLL_SEGRE_SUPPORT_SEPARATE_REVIEW.md)：元十二源の全被約支持はL6＋3×3rank一Segre五。全因子は一次×四次又は二次×二次×一次。旧四族whole-coverは棄却したままで、正しい全zero/unit補集合を保持 |
| 3 twistedcubiccone、rank0 | **全K又はØ** | [T38全源認識](CUBIC_CONE_ZERO_T38_WHOLE_RECOGNITION.md)：L3半単純の全W六を四actual shearで真completedT38へ、全L3非零上rank六kernelbundleと閉性でnil全fiberもT38。Q2全W三はH因子、base零はrank二以下 |
| 3 smoothscroll、rank2 | canonical dq2-root Ø | [A25 whole proof](CUBIC_SCROLL_A25_WHOLE_NONDOMINANCE.md)＋[erratum](CUBIC_SCROLL_A25_COORDINATE_AND_POLE_ERRATUM.md)：rank一/二、node nil/nil row pairとcusp固定点preserving pair、全ri・四fullP2 fibers。固定S像≤6から全≤24だけ。全pairのactual tensor正常化／親が残る |
| 3 twistedcubiccone、rank2 | base零rank≤2 | [CONE_RANK_TWO_REVIEW](CUBIC_CONE_RANK_TWO_SEPARATE_REVIEW.md)：全27、三Sym²base、Asem/Anil/BnilのwholeλとWsem/Wnil二fixed、全十九壁、vertex十八。五族の全E未了 |
| 4 rank0 | **全Ø** | [NIL32_ROOT_REVIEW](CI22_NIL32_AND_DEGREE_FOUR_RANKZERO_ROOT_REVIEW.md)：残single-eigen J3+2全六源のJacobi二係数−2a0²/−4b0²からa0=b0=0。全actualF=−q2²(a1²x2/2+a1b1x3+b1²x4)、零補集合はrank二以下。全centre/Jordan/Kronecker前段も戻して次数四全体が空 |
| 4 rank2、normalCI22／非CI scroll-cone | Ø | [NORMAL_CI22_REVIEW](NORMAL_CI22_RANK_TWO_SEPARATE_REVIEW.md)、rank4quadric Weil-class全排除 |
| 4 rank2、nonnormalCI22 | **全Ø** | [ALL_CI22_EXCLUSION](INTEGRAL_CI22_RANKTWO_IRREDUCIBLE_QUINTIC_EXCLUSION.md)：regular22_1はF=0又はH+2Q、nil32は全九inversechartと元compatibilityを戻しsmoothQ*C3又はH+2Q又はF0、M1は元Jacobiideal内H因子、node/cuspは全rank四空。全10/12/10/14/16源と全零補集合を保持。distinct-Q可約Dは§4へ戻り、ここに別のE穴を残さない |
| 4 rank2、projectedVeronese | **全Ø** | [SL2全支持の独立review](PROJECTED_VERONESE_TANGENT_SL2_SEPARATE_REVIEW.md)：全十五源Sym2+Sym4+Sym6、全Poisson cone次元二以下から全点は三Borel eigenlineへ。rank二／rank零空性／実Pf物理weight十二によるH因子で尽くす。これと全centre/CI22から次数四全rankが空 |
| 5 g0、rank0 | cusp／cone等Ø | [GENUS_ZERO_REVIEW](RANK_ZERO_QUINTIC_SURFACE_GENUS_ZERO_REVIEW.md)：normalization113/122、plane導手CP=2F等、CP選択＋全Grassmann投影源、固定D residue/Poisson一意。像22/24だけ |
| 5 g1、rank0 | **全K又はØ** | [不安定全源T40独立review](UNSTABLE_ELLIPTIC_RANKZERO_T40_SEPARATE_REVIEW.md)：任意E=A2⊕F(2,3)と全五sections/frameをstableへの逆拡大、相対MCM導手basechange、conductor係数descent＋Bottの一意ambient liftで実Poisson曲線へ。stableもunstableもT40、canonical/零summand/cusp等はØ |
| 5 g0、rank2 | cusp/oneprime node/plane-conductor不可能型Ø | [RANK_TWO_GENUS_ZERO_REVIEW](RANK_TWO_QUINTIC_SURFACE_GENUS_ZERO_SEPARATE_REVIEW.md)：残113/023 CP=Q、014 CP=Q+F等、全split二ΓL、全projection Schubert源13/9/3と全marks。像24/23/22だけ |
| 5 g1、rank2 elliptic | Ø | [ELLIPTIC_REVIEW](RANK_TWO_QUINTIC_SURFACE_GENUS_ONE_ELLIPTIC_SEPARATE_REVIEW.md)：全零summand／Pic0twist／node/cusp |
| 5 g1、rank2 canonicalindex2 | dependentmarks／cusp Ø | [INDEX_TWO_REVIEW](RANK_TWO_QUINTIC_CANONICAL_INDEX_TWO_SEPARATE_REVIEW.md)：全rank3q十一＋lowerframe三＋AB−u1C=q²の固定十一parameter、三monomial cases・α非零・LA/LB。真C*fiberで≤24。nodal全源の親未認識 |
| 6 rank0 | Ø | [SEXTIC_CONDUCTOR_REVIEW](RANK_ZERO_SEXTIC_CONDUCTOR_SEPARATE_REVIEW.md)：全normal/nonnormal S、五正規化scroll |
| 6 rank2 | 113/122/014、005非plane等Ø | [SEXTIC_RANK_TWO_REVIEW](RANK_TWO_SEXTIC_SURFACE_SEPARATE_REVIEW.md)：023全split-conductor finiteprojection proper locus、Gr(5,8)内dim≤14＋全frame、固定D closed-residue一意から≤24。proper locusの全components認識は別 |

degree4 rank0 nil32六源とrank2 nil32十二源は別入力で、各全源証明を用いている。両restriction rankの全integral CI22は既約Dでは空。
smooth nonplane rank零も全T40又は空。次数三rank零・次数四全rank・次数五genus一rank零は既全点認識済みであり、G上界だけの旧表を残余へ戻さない。
smooth cubic scroll rank零の旧四family whole-coverは独立反例で棄却のまま、正しいV2のL6＋全Segre五だけを用いる。
cut1044後のgenus零T40／cubiccone rank二T39／HF八壁／tangent m二等の進行中証明はここに混ぜない。

## 6. rank零plane：G完了が覆わない全first-symbol類

共通全源は[WHOLE_SUPPORT_NONDOMINANCE](PLANE_RANK_ZERO_WHOLE_SUPPORT_NONDOMINANCE.md)の35/30/9、
全multiplicityのabelian conormal。旧m二で除いたleading零parameterを全mへ戻す。

| 全first-symbol類 | K／Ø部分 | **Bに残る全源** |
|---|---|---|
| co-Higgs k1、rankM≥2 | M(p)=0 highmの全κopen/零はØ | [HIGH_M_REVIEW](PLANE_HIGH_M_RANK_AT_LEAST_TWO_SEPARATE_REVIEW.md)：rankM3 τ零sem/nil、rankM2 m2 kernel/nonkernelの全八family/walls、M(p)非零highm四源NDSS0/NDN0/BSS0(τ)/BN0(τ) |
| globalV、normalss、m二 | Cartan全periodはK=MPS | [M2_COMPLETION_REVIEW](PLANE_RANK_ZERO_M2_COMPLETION_SEPARATE_REVIEW.md)：mixed/regularnil/commonkernel等の二parameter四源・三parameter二源、全23plane/232curve/955root＋別6root等。各B層の元Poisson subsetの認識が必要 |
| normalnil b≠0 | 指定factor枝Ø | [ALL_MULTIPLICITIES_REVIEW](PLANE_NILPOTENT_AFFINE_ALL_MULTIPLICITIES_SEPARATE_REVIEW.md)：J3一orbit＋rank1全ρ/共鳴、highmρ−1/3 row零も含む。全像≤22だけ |
| lowparameter全21 | 一部factor／不整合Ø | [LOW_PARAMETER_CLOSURE](PLANE_RANK_ZERO_GLOBAL_VECTOR_LOW_PARAMETER_CLOSURE.md)：全一parameter線／onecut／compatibility残四・三linear枝／Tsquareの≤24。全actual親は別 |
| highm independent/common-image | K零二型はrank≤2でØ | [COMMONIMAGE_REVIEW](PLANE_GLOBAL_COMMONIMAGE_SEPARATE_REVIEW.md)：K非零二whole orbitと全higherは≤20だけ |
| dependent T=tV | constant-factor QはE8、残八Ore modulesはØ | [DEPENDENT_ORE_REVIEW](PLANE_DEPENDENT_ORE_SEPARATE_REVIEW.md)。ただし最初のH二源（normalss T0/K−5zw、normalnil T0/K5z²/2）は全像20/21のB。N三源のØ/E8をHへ外挿しない |
| traceless零、T≠0 | 指定F2因子枝Ø | [SCALAR_REVIEW](PLANE_RANK_ZERO_SCALAR_SEPARATE_REVIEW.md)：五Jordan×K ss/nil/zeroの十五源、全regularSS比・三δ0線・全根。≤23だけ |
| T零、K≠0 | Ø | 全kernel九・実shear/Eulerで全rank≤2 |
| weight零symbol全零 | anchor零／共通零はØ | [ZERO_SYMBOL_REVIEW](PLANE_RANK_ZERO_ZERO_SYMBOL_SEPARATE_REVIEW.md)：全39higher→A=O(−1)³ transitive Lie algebroid＋CE二cocycles九。aff1 kernel／P3 extension／全global cocyclesの≤24だけ。全orbit/algebroid正常形や親は未決定 |

上記finite wall数は重複込みcover record数で、成分数ではない。
次は保存された**全linear kernel＋元二次Jacobi＋全D**を使う。
Gに不要として解かなかったhigher compatibilityを、Eで暗黙に零又は自由としない。

## 7. rank二plane：Rとmを尽くした後の弱上界源

Rはintegral smooth/node/cusp、transverse/tangent conic+line、triangle/concurrent、2L+M、3L。
integral R全mは[IRREDUCIBLE_CUBIC_REVIEW](PLANE_IRREDUCIBLE_CUBIC_RANK_TWO_REVIEW.md)によりØ。m五以上もØ。

| 全枝 | K／Ø部分 | **B/Pに残るwhole source** |
|---|---|---|
| m二SF、全R | **全K又はØ** | reducedRはMPS、nonredRはØ又はE8。この枝へ未帰属源を足さない |
| m二double、transverse conic | Ø | 全irregular/moving/constantroot |
| 同tangent conic | h³全E8、一源bad-columnØ | [DOUBLE_TANGENT_COMPLETION](PLANE_M2_DOUBLE_TANGENT_COMPLETION.md)：h=f,Z0の固定源（K1八/K2零）、moving root全t≠0と全七壁は≤23のB。h³のV2恒等式を他源へ移さない |
| 同triangle | HF全Casimir源と大半H21はØ | [TRIANGLE_DOUBLE_REVIEW](TRIANGLE_M2_DOUBLE_SEPARATE_REVIEW.md)：残H21全十walls、K1七/八・K2零/一、真fiber九から≤23。全E未認識 |
| 同concurrent | 指定first-profile空枝Ø | [CONCURRENT_DOUBLE_REVIEW](PLANE_M2_CONCURRENT_DOUBLE_SEPARATE_REVIEW.md)：全17source/全QQ(α) finite pencil、全higher必要核、≤23。全17をMPSとしない |
| 同nonred moving／R3Lconstant | Ø／K | 全movingØ、R3L全八whole源はFULLP2九OreからE8 |
| 同R2Lconstant | 23fixed／大半wholewallsØ、最後nil全pØ/E8 | [CONSTANT_DOUBLE_REVIEW](PLANE_DOUBLE_LINE_M2_CONSTANT_DOUBLE_SEPARATE_REVIEW.md)：**八HF壁** k0,p5/8；k−1/2,p2/5；k−2,p−1/2；k4,p7/4、各ε0/1。全kernel和8/8,10/9,10/9,9/9、像21..23だけ |
| m三nonconcurrent reducedR | **全Ø** | [全二残源shear/P2独立review](PLANE_REDUCED_M3_WHOLE_SHEAR_FACTOR_SEPARATE_REVIEW.md)：Mの全七K1とSF-C3-unusedの全八K1を六shear＋残一/二係数へ定数全逆、元九P2・350/Pfから全D proper因子。旧七空源及び全root-cover逆を戻す |
| m三concurrent | **全K又はØ** | K0全E8、K22既nil22long23、K4a全23軌道は**K=T37**。[実PGL全150正逆・126Pfの独立review](K4A_T37_ACTUAL_MAP_SEPARATE_REVIEW.md)とclaim1020の全actual曲線・一軌道橋から、親24成分はT37に等しい。nonOreはØ。未列挙330根候補も全strictQ theoremが覆い、このwhole枝のE義務は完了 |
| m三SF、nonredR | **全K又はØ** | [normal-row独立review](PLANE_M3_INVARIANT_HYPERPLANE_SEPARATE_REVIEW.md)：最後のcdep C=a²,g2全十核のC三行零とFULL九P2からz因子で空。全十九first源/全Oreを戻し、他は既E8又はØ |
| m三generic21 constant00 | **全K又はØ** | [constant全源集約](PLANE_NONREDUCED_M3_CONSTANT_ROOT_POINTWISE_COVERAGE.md)：C-exact全三＋八壁は全normal-row因子で空。十四first源・二十二module・全壁/逆を戻し、残constantvectorは全E8へ |
| 同movingdouble10 | **全Ø** | [全源独立review](PLANE_M3_INVARIANT_HYPERPLANE_SEPARATE_REVIEW.md)：旧十一源・全壁・全Oreの補集合を保持。旧弱二源の全normal-rowと唯一τ−8/3元−2p8²行をFULL九P2 particular込みで戻してz因子、全空 |
| 同movingsimple01/02 | **全K又はØ** | [全七HB／六HAの実逆](PLANE_MOVING_SIMPLE_TWO_WHOLE_ORBITS.md)：37源/22壁/十moduleの初期17全源はzero-row又はτ−11/3元p8²から因子で空。HB全七→既Ga9-distinct23、HA全六→新moving-simple2/15 smooth23。m零閉包を保持 |
| m三generictriple、nonredR | **全K又はØ** | [WHOLE_POINTWISE_COVERAGE](PLANE_NONREDUCED_M3_GENERIC_TRIPLE_POINTWISE_COVERAGE.md)：90main/110wallの全200routeと全special/literal逆を戻し、八essential-three成分＋新二smooth23へ全点被覆。二Y wallも全九P1/全inhomogeneousP2から元p8² rowと全actualD因子で空。direct62だけに未帰属を残さない |
| m四、全R／全partition | **全E8又はØ** | [最後の六全shear独立review](PLANE_M4_SIX_WHOLE_SEPARATE_REVIEW.md)：reduced RのSF三＋211三は全六actual shearの定数逆でP1零、全九P2へのglobal二次補正後に全actualDがproper因子。旧nonred／rank一／22／4／31も戻して全位数四を被覆 |
| m四、全image-rank一／22／4 | **全K=E8**、SFならØ | [REPEATED_GEOMETRY](PLANE_M4_REPEATED_LEADING_GEOMETRY.md)：generic22/4はrank一へ必ず還元される。D=L g4+g5の全源を八essential-three成分が覆う。SFだけはordinary四重点のτ矛盾でØ |

corrected90のdirect62をG段階で自動的にKと読まず、全200routeの後続E集約と二Y全源因子証明を接続した。
m四もGの次元上界だけでKとせず、最後の六源に全shear/P2因子を接続して初めて全点被覆とした。
source名・全係数・全unit・normalizer逆は各receiptのproducer packetへexactに固定されている。

## 8. 次のEを有限問題へ分ける順序

これは計算指示又は新claimではない。既存の全source入口を失わず、何を証明すればその入口が閉じるかを列記する。
既存計画の[E_multiple_hyperplane／E_reduced_nonnormal](PLAN.md)に以下を枝分けするのが自然である。
別sourceへの全包含を証明した場合は、重複する義務を一つにまとめてよい。

| 優先lane | 再利用する有限化 | このlaneを閉じる最弱の結果 |
|---|---|---|
| E-ORBIT | nilshort r2/3、Ga9linear-w、stable二、plane1+1+3三 | whole sourceの既知親への実confluence、又は一般actual軌道の成分性／既新成分との全相違。単一pointの陰性を追加し続けない。K4a親はclaim1021で既T37へ認識済みなので再開しない |
| E-FINITE-SURFACE | rank二cubiccone五module、rank二quadricconeの保存済み全module | 元全Jacobiのreduced supportを有限whole正常形へ送り、全actualD/rank四条件を戻してØ又はKを証明する。次数四全rank／次数三rank零は既Ø/Kなので再開しない |
| E-FINITE-PLANE2 | m二tangent固定／moving壁、triangle H21十壁、concurrent十七源、HF八壁 | 全higher必要核と非斉次particularを同時に戻し、全parameter曲線又は全proper factorを出す。壁外のG上界だけから壁内membershipを延長しない |
| E-PLANE0-WALLS | rankM≥2、nil b≠0全ρ、non-Cartan globalV、lowparameter二十一、common-image K非零二、dependent H二、scalar十五 | 全保存sourceを既orbit／新componentへ有限化する。単なる各kernel dimensionや全β properをE結論に置換しない |
| E-GEOMETRIC-FAMILIES | A25 root pairs＋全ri、degree5 genus零projection、index二node、degree6 S23 split-conductor proper locus | 全parameter／marksの既知成分への包含を幾何で示すか、proper locus／bundle familyの有限strataに分ける。現在のGrassmann／bundle上界には、そのreduced componentsの全分類が含まれていない |
| E-LOW-MODULAR | rank一coutside、A0/c非零、A=c0最高部六閉包の全低次延長、横断[J,Q]=0 | 全constant／linear／quadratic sourceと元Jacobiを保持して空性又は親を示す。固定β評価の核から全orbitへの橋を省かない |
| E-ZERO-SYMBOL | plane Π0零の全Lie algebroid／CE二cocycles | 全actual algebroidとextension classの正常形、又は全tensor sourceへのknown-parent map。九cocycle係数とP3 extensionを局所anchorだけで捨てない |

先に有限fixed sourceやwhole一parameter族を認識することは合理的だが、最後の幾何family／low-modular／zero-symbolを未列挙のまま「E残余は数源」と縮めない。
E全完了の受理条件は、P_redの全分岐をK又はØに送る全点包含と、新component候補を採録する場合の成分性・全既存相違である。
固定pointの全D判定、有限jet、接空間、GBの支持診断だけではこの条件を満たさない。

## 9. G93入口の全34採録recordとのcrosswalk

番号は上記34-record snapshotの保存順であり、claim数やsource数の再定義ではない。
これにより、表§2–7がG閉鎖木のどの全sourceをEへ送り戻したかを監査できる。

| record番号／claimの略記 | 本稿の入口 | 全点で使える結論 |
|---|---|---|
| 1 NORMAL_CLASSIFICATION | §2正規D | K=T39 |
| 2 R4_QUADRIC_CUBIC_FULL_SUPPORT | §2 Q+C | K |
| 3 R4_THREE_HYPERPLANES_QUADRIC_FULL_SUPPORT | §2三H+Q | K |
| 4 G_SIMPLE_H_QUARTIC_ALL_GENERIC_EXCLUSION | §4 Hquartic | K/Ø、残Ga9linear-w B |
| 5 G_TWO_HYPERPLANES_CUBIC_ALL_GENERIC_EXCLUSION | §4二H+cubic | 後続全非plane K、plane三B |
| 6 G_H_TWO_QUADRICS_ALL_GENERIC_EXCLUSION | §4 H+二Q | plane等K、jointSSrank一／pure-nilのB |
| 7 G_MODULAR_NIL_RANK2_NONDOMINANCE | §3 nil31/22 | Kとshort r2/3／stable二B |
| 8 G_MODULAR_RANK1_CONSTANT_NONDOMINANCE | §3 rank一・A0/c非零 | K/B混合 |
| 9 G_UNIMODULAR_LINEAR_CASIMIR_NONDOMINANCE | §3 R13全延長 | Ø/B |
| 10 G_UNIMODULAR_QUADRATIC_PENCIL_NONDOMINANCE | §3 R22全延長 | Ø/B |
| 11 G_UNIMODULAR_PULLBACK_NONDOMINANCE | §3 S23と六閉包境界 | K/Ø/B混合 |
| 12 G_NONDEGENERATE_DEGREE_THREE_RANK_ZERO_NONDOMINANCE | §5 cubic scroll／cone rank零 | 後続1034/1035により空／T38 |
| 13 G_INTEGRAL_DEGREE_FOUR_RANK_ZERO_NONDOMINANCE | §5 degree四rank零 | 後続1022で全Ø |
| 14 G_DEGREE_FIVE_RANK_ZERO_GENUS_ZERO_NONDOMINANCE | §5 degree五g0 rank零 | projection family B |
| 15 G_DEGREE_FIVE_RANK_ZERO_GENUS_ONE_RECOGNITION | §5 degree五g1 rank零 | 後続1042により全K=T40又はØ |
| 16 G_DEGREE_SIX_RANK_ZERO_QUINTIC_EMPTY | §5 degree六rank零 | Ø |
| 17 G_DEGREE_SIX_RANK_TWO_ALL_NONDOMINANCE | §5 degree六rank二 | 大部分Ø、S23 proper projection locus B |
| 18 G_DEGREE_FIVE_RANK_TWO_GENUS_ZERO_NONDOMINANCE | §5 degree五g0 rank二 | Ø/B |
| 19 G_DEGREE_FIVE_RANK_TWO_GENUS_ONE_ELLIPTIC_EXCLUSION | §5 degree五g1 elliptic rank二 | Ø |
| 20 G_RANK_TWO_QUINTIC_CANONICAL_INDEX_TWO_NODE_NONDOMINANCE | §5 index二node | B |
| 21 G_RANK_TWO_QUINTIC_CANONICAL_INDEX_TWO_CUSP_EMPTY | 同cusp | Ø |
| 22 G_SMOOTH_QUADRIC_POISSON_CONORMAL_ALL_RESTRICTION_RANKS_HYPERPLANE_FACTOR | §5 smoothquadric両rank | Ø |
| 23 G_QUADRIC_CONE_RANK_ZERO_NONDOMINANCE | §5 quadriccone rank零 | 後続1025で全Ø |
| 24 G_QUADRIC_CONE_RANK_TWO_ALL_ZERO_DIVISORS_NONDOMINANCE | §5同rank二 | Ø/B |
| 25 G_INTEGRAL_DEGREE_THREE_RANK_TWO_ALL_NONDOMINANCE | §5 A25／cubiccone rank二 | canonicalrootØ、残B |
| 26 G_INTEGRAL_DEGREE_FOUR_RANK_TWO_ALL_NONDOMINANCE | §5 degree四rank二 | 後続1041により全Ø |
| 27 G_PLANE_RANK_ZERO_ALL_SYMBOLS_ALL_MULTIPLICITIES_NONDOMINANCE | §6全表 | 一部MPS/E8 K、factorØ、多数B |
| 28 G_RANK_TWO_PLANE_IRREDUCIBLE_CUBIC_EMPTY | §7 integralR | Ø |
| 29 G_PLANE_RANK_TWO_M2_ALL_NONDOMINANCE | §7 m二全表 | SF K/Ø、double K/Ø/B |
| 30 G_PLANE_RANK_TWO_MULTIPLICITY_FOUR_ALL_ZERO_DIVISORS_NONDOMINANCE | §7 m四全表 | 後続1040で全E8又はØ |
| 31 G_PLANE_REDUCED_NONCONCURRENT_CUBIC_M3_ALL_NONDOMINANCE | §7 m三nonconcurrent | 後続1043で全Ø |
| 32 G_PLANE_CONCURRENT_M3_WHOLE_NONDOMINANCE | §7 m三concurrent | 後続claim1021により全K/Ø |
| 33 G_PLANE_RANK_TWO_NONREDUCED_ZERO_CUBIC_M3_SQUAREFREE_ALL_NONDOMINANCE | §7 nonredSF | 後続1037で全K/Ø |
| 34 G_PLANE_NONREDUCED_M3_GENERIC21_ALL_NONDOMINANCE | §7 generic21三flag | 後続1036/1038/1039で全K/Ø |
| snapshot後1018 corrected-triple G | §7 generictriple | 後続1029で全K/Ø、plus-sign旧coverは不使用 |
| snapshot後1020 K4a smooth24 parent | §7 concurrent | 親24成分と全境界への包含、次の1021でKへ確定 |
| snapshot後1021 K4a=T37、wholeconcurrent全点被覆 | §7 concurrent | K=T37／既nil22long23／E8又はØ、新成分零 |
| 1022 degree四rank零全空 | §5 degree四rank零 | Ø、rank零nil32六源を全点で閉じる |
| 1023–1024 triple二weak wall成分 | §7 generictriple | K=新smooth23二成分、全prior126及び相互と相異 |
| 1025 quadriccone rank零全空 | §5 quadriccone rank零 | Ø |
| 1026 regular22_1 rank二全known | §5 CI22 | 既約DはØ、可約D全rank四は既H+2Q |
| 1027 triple二Y壁全空 | §7 generictriple | Ø、元九P1/P2とp8²から実因子 |
| 1028 cone-cusp rank二全rank四空 | §5 CI22 | Ø、全projective Casimirと両零補集合 |
| 1029 correctedtriple全点 | §7 generictriple | K=八essential-three＋上の新二23成分、又はØ |
| 1030 singularM1 rank二全既約D空 | §5 CI22 | Ø、元350ideal内H因子。可約distinct-Qは§4へ |
| 1031 nil32 rank二全known | §5 CI22 | 全十二源・九chart、既MPS／H+2Q又はrank二以下 |
| 1032 全integralCI22 conductor既約D空 | §5 CI22両rank | 全Ø、全pencilと全R／vertical補集合を保持 |
| 1033 moving-simple2/15新23 | §7 moving01/02 | K=新129番smooth23、全W六・全prior128/相互比較 |
| 1034 smooth cubic scroll rank零全空 | §5 degree三rank零 | Ø、修正L6＋Segre五だけを採用 |
| 1035 cubiccone rank零全T38 | §5 degree三rank零 | K=T38又はØ、全nil fiberも保持 |
| 1036 constant00全E8 | §7 generic21 | K/Ø、Cexact三全壁因子 |
| 1037 nonredSF cdep g2空 | §7 nonredSF | Ø、全K1/FULL九P2 |
| 1038 moving10全空 | §7 generic21 | Ø、元全壁とcritical二乗式 |
| 1039 moving01/02全known | §7 generic21 | K=Ga9-distinct／新HA又はØ |
| 1040 m四全pointwise | §7 全m四 | K=E8又はØ、六whole shears |
| 1041 次数四全rank空 | §5 全degree四 | Ø、PV全SL2分解＋全centre/CI |
| 1042 unstable elliptic rank零全T40 | §5 degree五g1 | K=T40、実sections/導手/ambient曲線 |
| 1043 nonconcurrent reducedR/m三全空 | §7 m三 | Ø、全六shear/全九P2/actualD |
| 1044 rank二plane/m三全pointwise | §7 全m三 | K/Ø、全R・全root・全higherを戻す |

93-file manifestに含まれない後続proofを無断採用したという意味ではなく、末尾はcut1044までに別途正式受理された補完である。
現在進行中又はcut後のgenus零／rank二cubiccone／HF壁／tangent m二等は本稿のK/B判定へ自動反映していない。
旧sourceの歴史的上界やwithdrawn coverを上書きせず、各行の最終受理scopeだけを取り出した。
