# 現在地 — 分類に必要な情報だけ

2026-11-07。active queueは `PLAN.json`、選択規則は `RESEARCH_POLICY.md` と `PROOF_ROUTE.md`。
目的はPGL5商前のPoissonスキームの**被約全成分の尽くし**。全軌道と全非被約構造は別段階である。

## 証明状態と保存状態

個別に証明・独立確認済みの下界は **111成分（一般rank四105、一般rank二6）**。
`research/COMPONENTS.json` が唯一の加算台帳であり、全上界・全G・全E・分類主定理は未証明。
現在の主張台帳は931件。plan/check/gateの成功は数学の認証ではない。

本bundleのGit親は `6ab28a327b64543ba0a7bf8cd302bea9def517e6`（104/916、187頁）。
今回の日本語主稿と201頁PDFは111成分・931主張を統合し、最終2passと主要頁可視確認を完了した。
このbundle自身のcommit SHAとorigin/mainはGitで確認する。
全八平面quinticの個別source/cusp・全周囲成分性・相違・全必要curve分類をそれぞれ独立確認し、
全essential3既約D支持の八成分被覆まで閉じた。保存後に完成したcone/concurrent等は次checkpointへ回す。

## 主義務

| 義務 | 現在の結論と残余 |
|---|---|
| G 一般beta尽くし | 未完了。未知beta優越成分の次元は25/26/27、一般の実Dは既約被約非正規quinticに限る。rank零支持は非plane全枝とplaneの大部分を非優越へ還元。rank二支持の残枝を保持 |
| E_multiple_hyperplane | 未完了。全H+2Q(Q既約)、nil22top、nil31long等は被覆済み。nil22image全十一源の残りはstable line-conic/three-lines二源。short nil31 r2/3、他の明示源も未帰属 |
| E_reduced_nonnormal | 未完了。三相異H+既約Q、既約Q+既約C3全支持は被覆済み。1+1+3の平面三配置/非plane nil、1+2+2の平面二配置/lambda0純nil、Hquartic型6とGa9 c0 linear-w、他の既約非正規支持が残る |
| A 成分監査 | 現111件は全周囲成分性と相違を保持。新候補だけを監査し、下界を総数として固定しない |
| M 主稿 | 日本語LuaLaTeX。受理済み証明を理解順に統合し、未証明の分類主定理を完成扱いしない |

## Gの幾何学的入口と受理済み支持

`DOMINANT_COMPONENT_GEOMETRY.md` を使う。正常な実退化因子は既存原稿
`normal_actual_divisor_recognition_ja.tex` と全独立監査により全点T39へ認識済み。
可約又は非被約な全実D型もG一般元からは排除済みであり、E全点被覆とは区別する。

rank零非plane支持は全てT40又は像次元24以下へ尽くした。
次数三scroll/cone、全integral CI22、projected Veronese、四projected degree四scroll、全coneを含む。
degree六はD正規化のsectional genus零・全五minimaldegree scroll・全導手node/cuspの証明で空。
degree五はgenus零の全projection上界、genus一のcanonical枝排除と全elliptic modelからT40又は24以下。
非normal支持、nef境界、non-split導手とcodimension二の補集合を保持した全証明であり、smooth/node一般値だけの結果ではない。
入口は `RANK_ZERO_SEXTIC_CONDUCTOR_EXCLUSION.md`、
`RANK_ZERO_QUINTIC_SURFACE_GENUS_ONE_RECOGNITION.md` と各別作者review。

**smooth quadricはrank零・二の全支持が空。**
`SMOOTH_QUADRIC_POISSON_CONORMAL_HYPERPLANE_FACTOR.md` と独立
`SMOOTH_QUADRIC_CONORMAL_FACTOR_SEPARATE_REVIEW.md` が全R可約/非被約、全vertical/square型を含む。
正しいEuler twistのnormal contractionsを同時global gaugeで消し、intrinsic/raw bracketの全補正からI²を得る。
全h-rowの次数とrank四gradientがspan hyperplaneのPoisson不変性、従ってactual h因子を強制する。
quadric coneではradicalが残り、この補題を流用しない。

rank二支持では全degree六とdegree五正規化genus零が空又は像24以下。
degree五genus一のelliptic全nef/GG modelは二つのample classの留数行列det10で空。
**canonical index二正規化は残る。**
`smooth/nonplane clean` の全枝はnode仮定を除いてT40又はprojected Veronese像22以下へ還元済み。
特異/verticalの残り、degree三・四のrank二全源、quadric coneは別義務である。smooth cubic scrollのrank二integral reduced Rは、normal quadratic全root/verticalと正規束c2八対九の幾何矛盾で空。全21源とrank二openの二回の弱いdimension判定は未閉で、R可約非被約/coneを続ける。

## rank零planeの現在の全分岐

全multiplicityでconormal Lie bracketは可換。rank零hypersurfaceが被約Pfに反することと
local modular flowの零曲面へのtangencyから従う。globalEuler liftをstabilizerと誤認しない。
これにより全weight零35・weight一30・weight二9源と全Lax条件を高重複度にも使える。

- **全m二源はMPS40又は像24以下。** coHiggsの全common saturation、Cartan認識、全21低parameterと残二/三parameter、全23plane/232curve/955root出現、ssnil六点を独立Fractionで確認した。
  `PLANE_RANK_ZERO_MULTIPLICITY_TWO_NONDOMINANCE.md` と `PLANE_RANK_ZERO_M2_COMPLETION_SEPARATE_REVIEW.md`。
- k一coefficientrank二以上は全mで24以下。全normalizer/全wall、rank三tau0、rank二M(p)零/非零とkappa-openの元Jacobi矛盾を保持。
  `PLANE_HIGH_M_RANK_AT_LEAST_TWO_SEPARATE_REVIEW.md`。
- normalnil b非零globalVは全mで22以下。唯一の高m共鳴row零も旧全kernelを戻して確認済み。
  `PLANE_NILPOTENT_AFFINE_ALL_MULTIPLICITIES_SEPARATE_REVIEW.md`。
- traceless零scalar源はT非零の全15正常形が23以下、T零K非零は全rank二以下。
  K零も含むPi0零の全39源はtransitive Lie algebroid A=O_P1(-1)^3、aff1 kernelのP3 extension、全CE二cocycles七以下から像24以下。
  `PLANE_RANK_ZERO_SCALAR_SEPARATE_REVIEW.md`、`PLANE_RANK_ZERO_ZERO_TRACE_BINARY_SEPARATE_REVIEW.md`、`PLANE_RANK_ZERO_ZERO_SYMBOL_SEPARATE_REVIEW.md`。
- 高m globalVの線形独立common-image四wholeorbitはK非零が20以下、K零は全18係数で共通partial_aとなりrank二以下。
  `PLANE_GLOBAL_COMMONIMAGE_SEPARATE_REVIEW.md`。

**残る高m dependent T=tV源は、全21Jordan族から五固定＋一wholeparameter線へ必要還元済み。**
全C/D極端normal Jacobi係数、actual hyperplane因子とquotient wedge lemmaで全例外を保持した。
`PLANE_RANK_ZERO_DEPENDENT_GLOBAL_VECTOR_HYPERPLANE_REDUCTION.md` と独立review。
六源の有限operatorで三源は20/21以下又は空となったが、残るK零三源の全G上界は未登録。
巨大rank二部分を含む線形上界33〜35を成功扱いせず、Ore suspension/元Jacobi/実essential変数数の構造へ移った。

## rank二planeの残余

零cubic Rが既約被約ならsmooth/node/cusp全型で空。
可約被約Rではm四squarefree leadingのgeneric j全nonconcurrent源をscalarへ還元し24/23/23以下。
特別jも全stable pencilを保持し、nonconcurrent三型を同様に処理した。
leading rank一squarefreeはD=L*g4+g5のordinary四重点mu9と全Saito源tau12の矛盾で全j/concurrentも空。
**concurrent特別j rank二pencil、leading重根、m二三、非被約Rは未解決。**

## 本質三変数の既約quintic

全strict canonical lift/Pfaffian CasimirからSaito(2,2) cotangent＋全六mag源へ還元済み。
全Fで磁気shear六行列が可逆、finite analytic typesと有限curveAutから有限個のactual23orbitへ還元した。
具体C45、F=x5+x²y³+y4zは新しい104番目の23全成分として独立確認済み。
true48cut/D101/M48、一normalの非零二次余核、全prior103との相違が証明書にある。

全八candidateのparam/implicit/Saito/free/fullGL2+mag源、他七のD102smooth又はD101非零二次はproducerで完了。
rootの全七ambient独立Fractionは5.542秒でPASS、別作者の全八explicit source/cusp認識も0.202秒でPASS。新七つは個別成分として採録済み。
全必要curve分類のpureproofも別receiptでPASS。全zero/basepoint/tangent/Cremonaexceptionalを保持し、全essential3既約Dは八成分閉和へ尽くされた。現在の直接分類証明ではSaito1971局所判定/未読Namba分類への必須依存も除いた。
入口 `ESSENTIAL_THREE_QUINTIC_EIGHT_COMPONENTS.md` と
`ESSENTIAL_THREE_QUINTIC_EIGHT_AMBIENT_SEPARATE_REVIEW.md`。現下界111/rank四105、全上界は未証明。

## Eの受理済み認識と次の数学

全H+2Q、三相異H+既約Q、既約Q+既約C3は全配置/全cofactor/全六magを保持して被覆済み。
1+1+3は十三平面配置のうち九新23成分と一toric衝突被覆、残三配置と非plane nilを保持。
1+2+2は七free配置の多くを既知族又は二新23成分へ被覆し、残二平面配置とlambda0純nilを保持。
詳細は主稿と `REMAINING_DEGENERACY_DIVISOR_SUPPORT.md`。G非優越だけの分岐をE帰属に置き換えない。

Hquarticのregular nil4、nil31long、nil22長translation、Ga9非零cと二根源の帰属/新成分は完了。
Ga9 c0 power1/power3の全源はA3+A2平面quarticと共通接Hの既知24成分へ。
**Ga9 c0 linear-wと型6 nil22imageの平面free quinticは残る。**
linear-wの全六normal二次支持はplane＋line/height四であり、原点根基や実normal曲線を得たとはしない。
marked A3+A2親族の全不可能性も保持し、三次以上の自動反復は停止。
型6はmarked quarticのMilnor/Tjurina/contactを使うwholecurve幾何へ移った。D5+tangent全markedpair＋六magは一actual23orbitへ還元済み、fullD101/true48cutの二次余核は零なので成分性未証明。高次へ自動反復せず、wholeparent又は他markedbranchを処理する。

nil22image全十一源のうち、triple二源は新23、double_B/Cと他の多重源は既知成分へ。
quadratic_A全七源はactual symplectic曲線と全period SNC親からMPSへ。
同じlambda7 anchorはsmooth唯一24、旧48cutの被約支持は一次元なので旧零次元目標は偽。
double_A全六magも、既知cuspidal24全源のa=1/4+epsilon正則族へ固定actualPGL/全mag可逆mapで入る。
**残るstable二源の全source閉包は23次元、個々のgeneric orbitは22/21。**
T01全閉包はorbit-rank20以下で両者を含まず、23次元既約D親族も非被約D閉条件と同次元からwhole stable源を覆えない。
この否定を全E帰属又は新成分性の証明とは呼ばない。

## 依存・再現・中止経路

現在の加算監査は `KNOWN_COMPONENT_AUDIT.md`、全結果は `research/CLAIMS.json`。
元134式/S61と全350Jacobi/150係数への同値・inverse・excluded chartsを維持する。
分母の補集合、generic supportと全point、計算された必要条件と全idealの同値を区別する。
旧source ZIPとimmutable producer/report/certificateは変更せず、compatible runner追記は旧FunctionDef ASTを保持し、旧実行bytesと契約hashをsnapshotへ保存する。

直接GBは五回未完了で停止。可逆交換子・cross・多項式左逆・平面五点・有限性だけでは第六点を排除出来ない反例を保存済み。
有限48/144忠実chartもunit未証明で、normal-only展開は密化した。
`NORMAL_DIRECTION_INCIDENCE.md` のF5+R表示や低次minorの再検算を、新たな尽くしとして繰り返さない。
B18のH0/T2/正確nil指数、Betti/接触精緻化、別稿publication sprintは休止。
CM/LPT/MPS/Pymの外部定理と旧版の全証明書を毎checkpointで再証明・再実行したとはしない。

mainへの通常commit/pushと現在の依存完全bundleを節目で保存する。生成・検証・GitHub保存・Library保存は別状態。
Libraryのupload receiptは無い。checkpoint後もprimary数学を続け、G/Eと上界が閉じるまで完了扱いしない。
