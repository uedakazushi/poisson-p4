# 日本語原稿・PDFの最終同期 — v61-s61

2026年11月3日（Asia/Tokyo）。[初回同期](MANUSCRIPT_SYNC_20261103.md)の後に、
独立監査で受理された非平方 d=0 全層の証明を追加した。一般五点の尽くしは未完了で、
既知成分数の下界57、総数の上界未確定を維持する。

## 最終版に統合した範囲

[主原稿](../papers/s61_first_jet_reduction_ja.tex)には、初回で統合した v60 の8項目と
非平方 a=0 全層の証明を保持し、d=0 の全層除去を独立の数学的な節として加えた。

- 元135式からの QQ[m] Smith 変換、係数と変換の p 整性、非定数不変因子を反転しない有限表示。
- 重み四561行269列の全先頭次数が零で、各列に非零定数 pivot があることによる全射。
- 一つの固定一次ジェット・支持面の上の第二因子の全 P¹と、実際の旗端点。
- 真の旗 core r=−4/25 における V-flat な平方零閉商と有限 kernel、端点の飽和支持に対する重みの議論。
- P¹の固有性と支持面の既約16次元因子への固有像を用いた、十分一般の β 上の全層への移行。

受理された証明と監査は [全層の証明](../research/notes/S61_BINARY_D_ZERO_STRATUM_REMOVAL_20261103.md)、
[混合標数の橋](../research/notes/S61_BINARY_D_ZERO_BRIDGE_20261103.md)、
[独立監査](../research/notes/S61_BINARY_D_ZERO_BRIDGE_REVIEW_20261103.md)を参照。
有限部は [fresh builder の成功記録](S61_BINARY_D_ZERO_PENCIL_20261103_V2/report.json)に対応する。

native 出力の `DET 1` は、269個の非零定数 pivot を除去した後の規格化された
0×0残余行列の値である。元の規格化前の269行基底の行列式を測定して1と得たわけではない。
全射を与える根拠は「全269先頭次数が零、全定数 pivot が非零」という計算であり、
原稿もこの根拠を記した。

[現行評価資料](../research/notes/p4_poisson_status_strategy_ja_20261030_v1.tex)、
[原稿索引](../papers/SOURCE_INDEX.json)、[原稿入口](../papers/README.md)も同期した。
二層の除去だけから全非平方 marked 境界を尽くしたとは推論しない。
残る不変因子の Jacobi 由来の排除、平方 c=2a の全周囲の解析、固定二座標候補の
全アフィン長と一般五点の尽くしは引き続き未解決として明記した。
係数模型の反例や中央商の枝を、実 Jacobi 解や真の境界枝として本文へ追加してはいない。

## 実測した最終成果物

|成果物|ページ|バイト数|SHA-256|
|---|---:|---:|---|
|[主原稿PDF](../papers/s61_first_jet_reduction_ja.pdf)|110|1210983|`3c536464d867bc44df9af7b331539ba591be1899bfdedfda0d8e18d33f2972c9`|
|[評価・現行追補PDF](../research/notes/p4_poisson_status_strategy_ja_20261030_v1.pdf)|14|451636|`f56987b4563e1389cf77708ca61cadd8a79d64a8c2b2ca612d6c8cea698d9dfc`|

両PDFを `pdfinfo` で測定し、フォーム・JavaScriptがないことも確認した。
ソース・索引・最終ログのハッシュは [最終機械可読記録](MANUSCRIPT_SYNC_20261103_FINAL.json)に保存した。

## 初回記録の保存

初回の109ページ版と当時の14ページの評価資料、TeX、索引、README、ログ、元の同期報告は
[最初の厳密なスナップショット](MANUSCRIPT_SYNC_20261103/first-snapshot-v61-a0.tar.gz)に保存した。
[保存受領情報](MANUSCRIPT_SYNC_20261103/first-snapshot-receipt.json)とアーカイブ内の
`SNAPSHOT_SHA256.json` にハッシュがある。アーカイブ内の全ソース・PDFが、元の
`MANUSCRIPT_SYNC_20261103.json` の記録と一致することを検証した。
初回の `.md` と `.json` は変更していない。そこでのパスとハッシュは初回時点の状態を表し、
現在の上書きされた作業ファイルではなく、このスナップショットと対応する。

## 最終検証と再ビルド

既設環境の `/workspace/cloud-setup/activate.sh` を読み、主原稿・評価資料をそれぞれ
LuaLaTeX で三回コンパイルした。セットアップや依存導入は再実行していない。
再ビルド方法は初回同期記録と同じで、生成したPDFを追跡対象の既存パスへコピーした。

最終ログに TeX エラー、未定義参照・引用、Overfull box はない。
既存の小型数式フォントの代替警告は残る。
新節のテキスト抽出に加え、主原稿103ページの d=0 証明と射程の段落を画像で点検した。
最終ログは [主原稿](MANUSCRIPT_SYNC_20261103_FINAL/main-pass3.log)と
[評価資料](MANUSCRIPT_SYNC_20261103_FINAL/status-pass3.log)に保存した。
v61 索引の6主張IDが `research/CLAIMS.json` に存在し、全参照ノートが存在することを確認した。
数学計算の再実行を、この文書コンパイル検査の範囲へ加えてはいない。

最終 native log は末尾空白を含む実出力を保持するため `.log.gz` に可逆圧縮した。JSONには圧縮前後のSHA256を記録した。初回版の実logは厳密スナップショットarchiveに含める。
