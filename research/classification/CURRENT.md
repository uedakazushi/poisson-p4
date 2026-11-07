# 現在地 — 分類に必要な情報だけ

2026-11-07。active queueは `PLAN.json`、選択規則は `RESEARCH_POLICY.md` と `PROOF_ROUTE.md`。
目的はPGL5商前のPoissonスキームの被約全成分の尽くし。全軌道と全非被約構造は別段階である。

## 証明状態と保存状態

個別に証明・独立確認済みの下界は **124成分（一般rank四118、一般rank二6）**。
`research/COMPONENTS.json` が唯一の加算台帳であり、全上界・全G・全E・分類主定理は未証明。
現在の主張台帳は980件。plan/check/gateの成功は数学の認証ではない。
最近の公開checkpointは `8b46f29ddbb9f539abc7a026e492ee25bf8ad699`（124成分・967主張、225頁）。
その後の受理済み純証明とfresh QQ replayを同じ234頁主稿へ統合し、通常main pushで保存する。
最新commit・push・PDFの実際の状態はGit及びbuild receiptで確認する。

Cloud中断後、初期snapshotへ戻ったworkspaceをfetch/ff-onlyで公開checkpointまで復元した。
失われた非公開reportを実測値ごと作り直さず、新REPLAY/RECOVERY_V2 identityで区別した。
exact手証明/sourceのhash復元と新算術は `../CLOUD_RECOVERY_20261107.md` と各manifestを参照。
旧任意Singular実体は消失し、二legacy producer gateがその旧binary identityだけでblocked。
既設M2・QQ/Fractionは動作し、影響を受けないprimary数学を続ける。セットアップ全体を再実行しない。

## 主義務

| 義務 | 現在の結論と残余 |
|---|---|
| G 一般beta尽くし | 未完了。未知beta優越成分の次元は25/26/27、一般の実Dは既約被約非正規quinticだけ。全rank零、rank二degree四/五/六・quadric、plane m四を排除。残るrank二degree三とplane m二三を処理する |
| E_multiple_hyperplane | stable line-conic/three-lines二源、short nil31 r2/3等の明示残源を実変形又は新成分へ認識する |
| E_reduced_nonnormal | E6smoothflex/A5hyperflex/Ga9linear-wの未知親、1+1+3/1+2+2の残源、他の非正規支持を全点認識する |
| A 成分監査 | 現124件の全周囲成分性と相違は受理済み。新候補があれば監査し、下界を答えとして固定しない |
| M 主稿 | 日本語LuaLaTeX。同じ原稿へ受理済み証明を理解順に統合し、分類主定理を未証明のまま明示する |

## G：全支持から残余を選ぶ

入口は `DOMINANT_COMPONENT_GEOMETRY.md`。全normal実Dは既全点T39へ、全可約・非被約Dも
未知優越一般元から除外済み。これはE全点被覆とは別命題である。

全rank零曲面支持は、非normal、nef境界、non-split導手、全multiplicityを含め完了した。
非planeはT40又は像24以下、planeはMPS又は像24以下である。
planeのcoHiggs/Cartan、scalar/zero symbol、common-image、dependent全mの全Ore八moduleを保持する。
`PLANE_RANK_ZERO_WHOLE_SUPPORT_NONDOMINANCE.md` と全別作者reviewを使い、旧粗33–35を再計算しない。

rank二の全smoothquadricは全Rで空、quadricconeは全Hamiltonianpencil/共鳴/compatibilityを含め像24以下。
全degree四はnormal CI22空、non-normal全五源像24以下、nonCI quarticscroll/coneのWeil class矛盾、PV像22で完了。
degree六及びdegree五genus零/elliptic全枝も空又は非優越。
**degree五canonical index二も全node像24以下、全cusp空まで独立完了した。**
全nonterminal/全従属marksを保持し、weighted symmetric normalization代数から全node11parameterと真C* fiberを数える。
cuspは導手idealとdegree四CIがdegree五Sを含めない短い矛盾で空。
入口 `RANK_TWO_QUINTIC_CANONICAL_INDEX_TWO_SEPARATE_REVIEW.md` と短cusp note。

smooth cubic scrollの全Rは負directrix C0を含む必要条件を得た。
normal canonical q2-rootも全rank/全R/vertex境界を保ちHartogs＋Casimir/divergenceの3対2で全空。
残るprimitive root A25=2C0+5fはrank一又は同class pairの二型だけ。
`CUBIC_SCROLL_DIRECTRIX_SEPARATE_REVIEW.md`、`CUBIC_SCROLL_CANONICAL_ROOT_SEPARATE_REVIEW.md` を使う。
A25 scalar connection及びquadric-net/P2-bundleは未受理の次routeで、四full-P2特殊fiberを除かない。
**smooth cubic scrollのA25二型及びtwisted cubic coneは未閉。** 旧21源/linear kernel16を全G上界としない。

## rank二plane

全integral reduced Rはsmooth/node/cuspとも空。
全m四は全R/root partitionを保持してG非優越まで完了し、nonredR全点はessential三の八23成分へ被覆。
m五以上はactualDがbinaryquinticとなり既約性に反する。

**m二・generic squarefree quadratic・全reduced Rは全点MPS40へ被覆した。**
全六primitive degree flags00/01/02/03/11/12を保持。
03/12はroot-line class矛盾、11 splitはactualPf係数矛盾、conic11は唯一marked位置のactualSNC曲線又は空。
constant-rootのtriangle/横断conicは全higher sourceのactualweight曲線でSNCへ、tangent/concurrentは
projective Jacobiのcurl `(2κI+4DY-DX)C=0` とnilpotenceから全空。
SNCに会う全ambient componentのMPS定理を使い、既知点の所属を曲線へ逆向きに外挿しない。
入口 `PLANE_RANK_TWO_M2_SQUAREFREE_REDUCED_CUBIC_COVERAGE.md`、全33file manifest、各独立review。
normaldouble及びnonredRはこの全点認識の補集合である。

m三R=3L squarefreeは既全rootcover/actualD3像で空。
R=2L+M squarefreeは全十九first symbolsを独立に尽くし、十六全枝は全leading imageのdiscriminant零から空。
c-dependent g二は実c→c+hb fiberによりB≤9、全P1十/P2 injectiveと合わせ像24以下。
fresh19 QQ replay22.097秒と独立4.412秒を別identityで保存した。
残DIAGscalar/NIL_SPLITはOre全Q＋Wへ数学的に変更し、全四module producerは空又はessential三まで完了、
独立確認を待つ。この段階で全m三G/Eを完了扱いしない。
reduced reducible Rのm三、全repeated leading、m二double/nonredRは別の残義務。

## E：全点認識と未知親

全H+2Q、三相異H+既約Q、既約Q+既約C3は全配置/全source/全magで被覆済み。
全essential三変数の既約quinticも全必要curveリスト・source・ambient・相違を独立確認し、八23成分へ全点被覆。
marked μ六九新23成分、μ五四新23成分を採録し、全下界124を保つ。
D5cusptangent全源は旧A0、3A2cusptangent全源は旧DPM、全A3A2共通接Hはk=-1も含め既cuspidal24へ。

全global T39/T40のrank四評価fiberはconnected Sp4とglobal Ga4切断で既約、S61交点は閉A16 graphと一致。
E6smoothflex、A5smoothhyperflex、Ga9c0linear-wの同一旧pointを別Fractionで全150両方向/134/350/45確認し、
両全global親の外まで証明した。全sourceのdense23orbitに適用し、zero-invariant境界は除外しない。
未知親・次元・新成分性は未解決のまま保持する。

T37の固定Lagrangian B枝は有限surjectiveとdegree20 field wallから全sourceがgenericlift閉包と一致し、
全actual restricted divisorがLQ²。上記三source及びstable二源への**指定uniform normal-cone入口**を除外した。
全T37 componentや他の座標縮退からの親族性は排除していない。未復元nativefactorを新再計算としない。

nil22image残stable二源の全source閉包は23、generic orbitは22/21。
T01と既23既約D親のwhole包含は排除済みだが、全E又は新成分性は未証明。
1+1+3平面三配置/非plane nil、1+2+2 plane二配置/lambda0純nil、short r2/3も残る。
詳細は `REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` と主稿。

## 再現と保存

元134式/S61と全350Jacobi/150係数への同値・逆写像・全局所化補集合を維持する。
旧source/report/manifestは不変、compatible追記は旧FunctionDef AST及び実行bytes/全契約を凍結してからcurrent bindingを更新。
直接GB五回とnormal-only密化、旧minor/jet反復、B18 H0/T2/nil指数/別稿仕上げは休止する。
各新finite jobをgate後実行し、負結果を証明にしない。一次文献の適用scopeは各証明で確認する。
通常main commit/pushで保存後も有効なprimary数学を続ける。Library upload receiptは無い。
