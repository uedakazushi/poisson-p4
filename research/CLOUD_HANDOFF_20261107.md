# Cloud研究引継ぎ — 現行分類を継続する

2026-11-07。公開済みCloud、repository `/workspace/poisson-p4`、branch `main`。
移行先で継続中の研究であり、旧移行一時停止を研究停止命令として扱わない。

## 保存境界と復旧

最初のCloud開始SHAは `4c3f623fb698f9adc6077015b25e2325b11072af`、
11月7日の数学baselineは `a00a541b835b2bff18bafcba545ecbd4bae7175b`。
復旧区間の開始SHAは `8b46f29ddbb9f539abc7a026e492ee25bf8ad699`、直近公開保存は
`57f4b36c00071103e4288793e1e3470f2ccf8e2f`（124成分/980主張/234頁）。
旧8b46は124成分/967主張/225頁。現行公開57f主稿は234頁、今回126/991の新保存単位は244頁。Cloud中断後に初期snapshotへ戻ったworkspaceを
fetch/ff-onlyでここまで復旧した。復旧時の8b46/PDF hashはCLOUD_RECOVERYに保持し、現公開57fの234頁PDF hashは
`fd82f5ed1a82e72b86a8b1ebedbcca7a27d5057438ae6e226f6b3a7cc7018451`。

非公開部分の旧report実測・旧hashを捏造して復旧しない。手証明/入力/sourceは旧hash一致のものだけexact復元と呼ぶ。
新m三算術はNEW_REPLAY（22.097秒）、独立算術4.412秒、splitreview metadataはRECOVERY_V2の別identity。
詳しくは `CLOUD_RECOVERY_20261107.md` と各manifestを使う。
次の意味ある受理証明を通常main commit/pushで保護し、その後もprimary数学を続ける。
実際の最新push状態・原稿頁数はGit/build receiptで確認し、実行前に成功と記さない。

## 現行数学

下界は126（rank四120、rank二6）、主張991。E6smoothflex/A5hyperflexの全二法線mprimary支持から二新23全成分を独立受理した。
全rank零、rank二degree四五六・quadric、plane m四全RをGから排除済み。
全essential三変数既約quinticは八23成分、μ六十/μ五五markedquartic新成分も個別ambient/相違まで受理済み。
新たにm二SF/reducedR全点をMPS40へ、degree五canonicalindex二全node像24以下/cusp空、
scroll全RにC0必要/canonical q2root空まで別作者確認した。
E6/A5/Ga9linear-wのdense23orbitは全globalT39/T40外、全sourcezero境界は保持。
T37固定B曲線の全finite wall閉包/LQ²は指定uniform入口だけを排除し、全T37親は排除しない。

全G/E・上界は未証明。smoothcubicScroll A25全二型も像24以下へ完了、coneも全27source/五module/十九壁の独立確認で全G非優越が完了、残るGはrank二plane m二三の残型だけ。
m三nonredR/SFは十九有限化の十六空＋cdep像24以下＋残Ore四moduleの独立確認を受理、全G像24以下が完了。
Eはstable二源、Ga9未知親、short r2/3、1+1+3/1+2+2の明示残源等。
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
