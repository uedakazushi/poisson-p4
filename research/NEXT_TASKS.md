# 次の課題 — 全被約成分の尽くし

2026-11-07 / v64-classification-first。現行キューはclassification/PLAN.json。

```sh
.venv/bin/python scripts/classification.py check
.venv/bin/python scripts/classification.py next
```

**個別成分の下界103、一般rank四97、rank二6。上界・一般五点・全分類は未証明。**
ユーザーの11月7日の再開指示により実作業を継続する。保存や旧移行用一時停止を終了条件にしない。
G、E_multiple_hyperplane、E_reduced_nonnormalを並行し、意味のある進展は主原稿・引継ぎへ統合して通常main pushで保存する。

## G：既約被約非正規quinticだけが残る

未知beta優越成分は射影次元25〜27。全非被約因子、全可約因子の次数分割、正規因子は
既知成分被覆又は全source次元上界で未知一般元から排除済み。
全1+1+3ではB0非平面8/7源を含む。全H+quarticでは型6のSaito次元論、型7の全liftを含む。
このG排除とEの全点帰属を混同しない。入口はDOMINANT_COMPONENT_GEOMETRY.mdと
REMAINING_DEGENERACY_DIVISOR_SUPPORT.md。

既約非正規五次の曲面に対し、今回次を証明確認した。

- smooth、generic ordinary node、normal-cone交叉にvertical divisor無しなら、rank0はT40、rank2は唯一projected Veronese。
- 全PV tangent源は十五変数、全350Jacobi/10079 liftの独立確認で固定源次元≤2、全PGL射影像≤22。rank0PVは非零Poisson無し。
- normal曲面のrank2枝では解消κ≥0がconormal flat moduleに矛盾し、stable楕円scrollのrank2も空。
- hyperplane内degree≥3曲面は全rankで既約五次に不可能。smooth CI22のrank0全枝とrank2-node全vertical枝も空。
- integral非退化曲面のlinear-normality例外は唯一smooth PV（特異Severiの一次射程確認済み）。Euler全sourceによりq=0・degree≥5の全rank0枝は空。

normal degree五rank0はTironi/Schröer分類とEuler完全源から全T40へ認識済み。CI22の全半単純pencilも実D零又はlinear因子付きと証明済み。
degree三全二軌道はabelianconormal全源と完全dimension証明でG像≤24、regular CI22全Jordan型もG像≤24を確認済み。
all-degenerate CI22も全九型・全二十源で空と独立確認し、全integral CI22のG像≤24が完了。全smooth非plane rank零もclean/node/vertical無しにT40又は像≤23へ還元した。
全degree四rank零は全centre/extraclassを保持した証明でG像≤24、normaldegree六も全nef数値式とrelativeconicで空となった。degree五六の全coneも全Euler sourceのQ因子で空。
D正規化のsectional genus零→minimaldegree全scrollと全導手を独立確認し、非normalも含む全degree六rank零を空とした。degree五のgenus零（残余plane付き）も全projection像≤24。
degree五genus一もcanonical排除・全elliptic bundleのT40/像≤24認識を受理し、quadric cone全源もOreと全一変数例外処理で像≤23。これで全rank零非plane枝がGから排除された。
rank零plane位数二k1 rankM≥2は全像≤23、Cartan全higher源はMPS40へ。全degree六rank二とdegree五rank二genus零は像≤24、種数一ellipticは全nef/GG型を排除。
全重複度conormal可換、normalnil b非零全源≤22、全essential3既約D≤23も受理。m4squarefreeleading rank一はordinary四重点の矛盾で全j/concurrentも空。
次はrank零planeのglobalV残nonCartan/commuting nilnormal・位数三以上、degree五rank二canonical index二、他のrank二特異/vertical支持、rank二planeの零R可約/非被約を処理する。smooth非plane cleanは全横断型を処理済みで、vertical/特異Sを残す。degree3の全Jacobi二次rowspaceは実Dを含まないので、
この最初の線形判定からrank四不存在を主張しない。旧五回未決GBやjet/minor列へ自動復帰しない。

rank二planeの既約被約Rはnode/cuspも含め全排除済み。非planeの三重以上はsecant論で不可能、rank零first Lieの全二型からconormal abelianを得る。この線形必要条件でdegree三全源を縮小しており、旧30変数GB/冪探索を繰返さない。

## E：未帰属sourceの実変形又は新成分性

H+2Q（Q既約）、五相異なるH、三相異なるH+既約Q、既約Q+既約C3の全支持は被覆済み。
rank2の全支持も六成分とrank4閉包へ接続済み。これらを欠陥なく再分類しない。

H+既約quarticは全半単純24源、純nil型1/2/3/4/5/8/10を完了した。
型4の唯一六源はliteral P+eR+e²Tでregular semisimple型12へ移り全MPS。
型5は新22成分とT40、Ga9型9のc=1全十源は新23成分へ被覆した。
型7全七源は三新成分（二22と一23）、MPS及び既存既約五次23成分へ全被覆済み。Hquarticに残るのはnil22image型6とGa9 c0 linear-w一源である。
Ga9 distinct-c全源も実cubic曲線で新しい既約非正規quinticの23次元軌道成分へ被覆済み。Ga9c0 power1/power3全源も新平面Hquartic24成分へ被覆済み。残るGa9はc0 linear-w一源だけ。全六normalの二次初期支持はplaneとline（height4）で原点に閉じない。診断Macaulay両包含を保存し、高次jet反復を止めactualwhole変形又は支持幾何へ切り替える。

重複Hではnil31long-top六源、nil22top八源、nil31image half/third全源は被覆済み。
nil22imageの全十一源は独立606恒等式で有限化した。triple_A/Bは二つの新23成分、
quadratic_B/triple_C/double_Dは既知25成分への全parameter実曲線で被覆。
double_B全六源もCPLUS_EXT24へ実変形で被覆済み。double_C全六源も従属HHHQ24へ全parameter実曲線と全八源actualmapで被覆済み。残る三源はstable二源、double_A。quadratic_A全七源はactualsymplectic曲線で全MPSへ被覆済み。double_Aの既知候補はT36/T37/T38まで限定したが帰属未証明。
short r=2/3も未帰属で、三次零を高次反復せず幾何又は実変形を使う。

1+1+3には平面三配置と非平面nil源、1+2+2にはplane tan-simple/chordとlambda零のE残源がある。
全G次元上界だけではこれらを被覆したことにならない。既約非正規quinticの全点帰属も残る。

## 計算・保存の規律

非自明計算はactive義務への矢印、全環・全補集合、有限資源と証明書をjob契約にし、gate後に実行する。
同経路二回未決なら数学的障害を診断して変える。最弱の支持・次元・actual curveを優先し、
正確nil指数や全非被約環を分類の隠れた前提にしない。B18 H0/T2と旧minor farmは休止。
Mは同じ日本語主原稿へ統合する作業、F/Pは全穴が閉じるまでblockedのままとする。
