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
その後の通常main pushで `9468048f7a9daa051fa57054858aa8f07a1dd61d` を公開し、HEAD=origin/mainを確認した（126/1020/288頁、一般G完了）。
旧8b46は124成分/967主張/225頁、旧57f主稿は234頁。直前136aの126/991保存単位は244頁。今回126/1002主稿は256頁。Cloud中断後に初期snapshotへ戻ったworkspaceを
fetch/ff-onlyでここまで復旧した。復旧時の8b46/PDF hashはCLOUD_RECOVERYに保持し、直前136aの244頁PDF hashは
`1df0210441101fbd577fc6278cbaff21ff800c36a7ccb9aa3b4fcda7096fe9fa`。

直前256頁PDF SHAは `5e95c800fb86f37f0cd95b6be1067e4a97bdd7f5d15f2a8058a5f91962528112`。
直前288頁PDF SHAは `a6482f56471bd0e75ae35a01a234a81d00ef30b7949d372cca910415da65e5ba`、TeX SHAは `9155f735d319067d99650d84159979e2572d3b939013b4f184b59d16ab6d6905`。
直前295頁PDF SHAは `b223a24182b809cc7dc0a2de6b5e2247bf80ea02070a72ae22d893477f43505a`、TeX SHAは `e9269978ea4239205676b1395b2e90635bddc6f8edb2a60d567d31197025bdd2`。
今回の継続区間開始SHAは `57f8e29db085dd43df885f9085b34ae6ea9053fa`、HEAD=origin/mainを確認した。
現行129/1044・301頁PDF SHAは `921b112e2d90aca3327317332a14420b6d19d3a51af61142c64deb3cbfb66b6b`、
TeX SHAは `340bf6035664c3f031ef62215243f6110b0420a4a88192e7c28949715d097b97`。
公開保存の実SHAはGit HEAD/origin/mainで確認する。

非公開部分の旧report実測・旧hashを捏造して復旧しない。手証明/入力/sourceは旧hash一致のものだけexact復元と呼ぶ。
新m三算術はNEW_REPLAY（22.097秒）、独立算術4.412秒、splitreview metadataはRECOVERY_V2の別identity。
詳しくは `CLOUD_RECOVERY_20261107.md` と各manifestを使う。
次の意味ある受理証明を通常main commit/pushで保護し、その後もprimary数学を続ける。
実際の最新push状態・原稿頁数はGit/build receiptで確認し、実行前に成功と記さない。

## 現行数学

採録cutは129成分（rank四123、rank二6）・1044主張。
新HA moving-simple t=2/15の全W六を一つのsmooth23全成分へ認識し、全prior128との相違を独立確認した。
rank二planeのm三・m四は全R/root/normal/higher源を既知成分へ被覆した。
smooth cubic scroll rank零は修正L六＋Segre五の全源で既約quintic枝が空、cubic cone rank零は全T38。
degree四導手曲面は全rank零/二・全centre/spanで既約quintic枝が空。
不安定elliptic rank零は全sections/frame/conductorの実曲線でT40へ被覆した。
一般Gは完了、全E・総上界・分類主定理は未証明である。
現行PDFは129/1044・301頁、最終二回組版と全参照・主要頁可視確認済み。
通常Git保存はHEAD/origin/mainを確認する。全H+2Q、三相異H+既約Q、Q+C3、1+1+3非plane等は既被覆済み。
旧smoothscrollfour-coverは反例により撤回、旧bytesを保存し正しい二源だけを採用した。

次cut用のgenus零rank零T40、cubic cone rank二T38/T39、HF八壁T36/T38、tangent m二、
dependent-global-vector、triangle/concurrent及びquadric cone rank二等のpacketは、独立確認と台帳採録を区別して保存する。
現1044件へ未採録のpacketを既証明の全分類へ算入しない。
genus零の旧未コミットreview note一件は同名作成衝突でbyte消失した。旧hashを改変せず、
`classification/inputs/genus_zero_rankzero_review_source_loss_20261108.json` に欠損と旧三記録のexact bytesを明記した。
原証明と新別名の独立reviewは現存し、欠損旧noteを証明依存には使わない。

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

全rank二plane m三及びm四のpointwiseEを閉じた。全degree四導手もrank零/二で既約quintic枝は空、unstable elliptic rank零全源はT40へ。
