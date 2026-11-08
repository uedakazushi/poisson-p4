# Cloud研究引継ぎ — 現行分類を継続する

2026-11-08。公開済みCloud、repository `/workspace/poisson-p4`、branch `main`。
移行先で継続中の研究であり、旧移行一時停止を研究停止命令として扱わない。

## 保存境界と復旧

最初のCloud開始SHAは `4c3f623fb698f9adc6077015b25e2325b11072af`、
11月7日の数学baselineは `a00a541b835b2bff18bafcba545ecbd4bae7175b`。
11月8日の直前継続区間開始SHAは `136a2a3ccb459da913b5db2e274bb5389830b9a6`。
今回の継続区間開始SHAは `dab368089679508ae5c9016d66ba6accfe6d8fa9`、GitHub mainとの一致も確認した。
復旧区間の開始SHAは `8b46f29ddbb9f539abc7a026e492ee25bf8ad699`、直前の公開保存は
`dab368089679508ae5c9016d66ba6accfe6d8fa9`（126成分/1002主張/256頁）。
旧8b46は124成分/967主張/225頁、旧57f主稿は234頁。直前136aの126/991保存単位は244頁。今回126/1002主稿は256頁。Cloud中断後に初期snapshotへ戻ったworkspaceを
fetch/ff-onlyでここまで復旧した。復旧時の8b46/PDF hashはCLOUD_RECOVERYに保持し、直前136aの244頁PDF hashは
`1df0210441101fbd577fc6278cbaff21ff800c36a7ccb9aa3b4fcda7096fe9fa`。

直前256頁PDF SHAは `5e95c800fb86f37f0cd95b6be1067e4a97bdd7f5d15f2a8058a5f91962528112`。
現行288頁PDF SHAは `a6482f56471bd0e75ae35a01a234a81d00ef30b7949d372cca910415da65e5ba`、TeX SHAは `9155f735d319067d99650d84159979e2572d3b939013b4f184b59d16ab6d6905`。
公開保存の実SHAはGit HEAD/origin/mainで確認する。

非公開部分の旧report実測・旧hashを捏造して復旧しない。手証明/入力/sourceは旧hash一致のものだけexact復元と呼ぶ。
新m三算術はNEW_REPLAY（22.097秒）、独立算術4.412秒、splitreview metadataはRECOVERY_V2の別identity。
詳しくは `CLOUD_RECOVERY_20261107.md` と各manifestを使う。
次の意味ある受理証明を通常main commit/pushで保護し、その後もprimary数学を続ける。
実際の最新push状態・原稿頁数はGit/build receiptで確認し、実行前に成功と記さない。

## 現行数学

下界126（rank四120、rank二6）、主張1020。**一般β Gは完了**、全E・総上界・分類主定理は未証明。
直前公開主稿は1002cut/256頁。現行保存単位は1020cut/288頁、G完了とK4a完全有限切片まで統合済み。

Gは最終 `classification/GLOBAL_GENERIC_G_EXHAUSTION_CLOSURE_AUDIT.md` と93入力manifestまで独立受理し、完了した。
全actualDの因子分割、非正規導手曲面rank零/二・全degree1–6、plane全m/R/normalroot/全higherを尽くす。
exactC0 bridgeの未知優越成分25–27は既知閉和、global像24以下又はproperβ像に入れず、残余B=0。
一般係数体の全fiberは既知五点。新一様分母/単位ideal展開を実行したとはしない。
旧plus符号90wholecoverageは撤回のまま。correctedminus新90・valid14joint逆・18modules/34walls・newfixed原ideal因子を三phase独立受理。
二E1DX壁はconsistentW6/像23以下で、target実在とEは未知。全concurrent m三も全K0まで独立完了。
解決済みGを欠陥なく再開せず、全E/成分帰属と主稿へ作業を進める。

EではH+二相異既約quadricのplane全七配置を全磁気で既知成分へ被覆済み。
1+1+3の全非plane支持も、全11nil表とj21(h,x)/j21(0,h+x)/p0の三actualcurve及び全n6零補集合から、
全MPS40又は旧同22成分へ被覆した。型の残Eは正確にplane三配置全六磁気源だけ。
他の残Eは1+2+2 nonplane/lambda零、stable二源、Ga9linear-w未知親、shortnil r2/3等。
全essential三変数既約quinticの八23成分被覆、現126成分の全ambient成分性/相違は独立受理済み。
新候補があれば126を答えに固定せず監査する。

K22 whole一軌道は既nil22long23 denseorbitへ全十actual括弧/全逆で同定し、全E所属も完了。
K4aは同じD101/G48からnormalizedC*全IFTgraphを次数5へ強制し、元350の全次数10零を別Fractionで確認。
真の102切片と48実作用のétale積でsmooth unique24、整数actual曲線の全非零parameterは同一dense24軌道。
全K4a23軌道はこの24成分境界へ。既知親名/相違は未認識、新加算零。入口 `K4A_COMPLETE_TORUS_SLICE_SEPARATE_REVIEW.md`。

1020cut後のE資料は各自独立receipt付きで保存中：K4a実24親は全order二零点operatorでT36/T38を排除、T37又は別24が残る。
二E1DX弱壁は全W6を実unipotent逆で各一23軌道へ縮約するpure証明があり、親の認識は未了。
shortnil r2/3はmarkedSS2+c closureの新障害と、一回のSS12 spectralframe endpoint不整合を保存。
これらを1020台帳/PDFの追加定理として黙って数えず、次のroot読解・採録を待つ。
全G-bound sourceの全点E inventoryも並行監査中で、上記数源だけで全Eを尽くしたとはしない。

全詳細と次の依存入口は `classification/CURRENT.md`、`classification/PLAN.json`、`NEXT_TASKS.md`。

## 環境と実行

既設 `.venv`、`/workspace/cloud-setup/activate.sh` を再利用。
Macaulay2 Cloud1.24.11対元1.26.05のQQkernel/既定GB互換小例は確認済み。全setupをやり直さない。
文献接続/一次scope確認も既receiptを再利用し、未取得Saito原全文を取得済みとしない。
独立算術はSymPy1.14又はFraction。日本語LuaLaTeX、PDF可視、通常fetch/pushは動作済み。

**旧work内の最適化Singular実体は復帰snapshotに無かった。**
旧任意binary SHAを要求する二legacy gateはblocked、QQ証明の失敗ではない。
別binaryを旧SHAに見せず、必要なら既設M2又はQQでfresh identityの計算を実施する。
既存`legacy/`はtracked原ZIPの展開入力でありforce-addしない。必要memberだけexact hashで回復する。

root統合担当だけがGit/大域台帳を編集し、log_manuscript担当が同じ日本語主稿を編集する。
旧FunctionDef/producer/report/manifestは不変、compatible追記前に旧実行bytesと全契約を凍結する。
各finitejobは義務/矢印/全環/全補集合/有限resourceを記してgate後実行。checkerは数学を検証しない。
二回非閉鎖なら数学routeを変え、GB/minor/jet列を増やさない。B18等secondaryは休止。
通常main pushはユーザー許可済み。Library upload receiptは無い。
現在の数学作業にユーザーの追加許可を要する阻害は無く、残る数学的義務を継続する。
