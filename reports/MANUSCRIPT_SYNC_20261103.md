# 日本語原稿・PDFの同期 — v61-s61

2026年11月3日（Asia/Tokyo）。公開済み Cloud 環境の既設 LuaLaTeX を再利用した。
セットアップ・依存導入は行っていない。

## 統合した内容

[一次ジェット還元原稿](../papers/s61_first_jet_reduction_ja.tex)へ、保存済み v60 の8項目を統合した。
三次元基本サイクル `65[B]+5[Fβ]` と次数30/70、一般切断長2300と t=0 商2270、
任意の固定二座標での t² コロンによる零次元候補、191生成表示による t 可逆長の回収、
有限 marked 点の剰余次数10以上と仮想曲面の交叉数20以上を一つの数学的な節にまとめた。
指定平方旗の純三次端点因子・h=0因子、c=2a の全有限パラメータ多項式模型と射影端点も記載した。
既存の零対層の証明には独立レビューを接続した。

独立レビューで受理された [非平方 a=0 全層の証明](../research/notes/S61_BINARY_A_ZERO_LIFT_20261103.md)を別節に統合した。
元135式の中央支持、同じ incidence の一般重複度と特殊横断上界、混合標数の自由加群、
飽和支持と非正重み底上の有限性、正しい c=0 の射影端点、P² の半不変因子の交叉を区別している。
指定 β で微分階数15であることを、同じ incidence の一般 β 優越性と混同していない。
V=Z_(32003)、完成環の係数環 Z_p、p-flat 閉包 `(I:p^∞)` の役割を明記した。

[評価資料](../research/notes/p4_poisson_status_strategy_ja_20261030_v1.tex)には現行追補を加え、
2026年10月30日の第24版評価を歴史的な本文として保存した。
本文の当時の未実行課題を現在の未解決課題として読むことがないよう、追補に優先順位と現行入口を記した。
[原稿索引](../papers/SOURCE_INDEX.json)と[原稿入口](../papers/README.md)も同期した。

一般五点の尽くし、d=0 の全層除去、平方 c=2a の全有限値での除去、特殊ファイバー全分類、
新規性の網羅的調査は本同期の完了項目に含めない。既知成分数の下界は57のままである。

## 実測した成果物

|成果物|ページ|バイト数|SHA-256|
|---|---:|---:|---|
|[一次ジェット還元PDF](../papers/s61_first_jet_reduction_ja.pdf)|109|1202141|`f172494288bc11a8c689441890941aad338e29fabf0e836acb9615da82ed862e`|
|[現状評価・現行追補PDF](../research/notes/p4_poisson_status_strategy_ja_20261030_v1.pdf)|14|449920|`da3d9bda19b5cae063185580eeb7dc00cd26aafd521dec77a8b6f1c454835cd2`|

旧104ページ・12ページという値は引き継がず、最終PDFの `pdfinfo` で測定した。
両PDFにフォーム・JavaScriptはない。ソース・索引・最終コンパイルログのSHA-256は
[機械可読記録](MANUSCRIPT_SYNC_20261103.json)に保存した。

## 再ビルドと検証

リポジトリ直下で実行する。

```sh
source /workspace/cloud-setup/activate.sh
mkdir -p build/papers
for n in 1 2 3; do
  lualatex -interaction=nonstopmode -halt-on-error \
    -output-directory=build/papers papers/s61_first_jet_reduction_ja.tex
done
bash scripts/build_notes.sh
```

実際に両文書を三回ずつコンパイルし、最終PDFを追跡対象の既存パスへコピーした。
最終ログに TeX エラー、未定義参照・引用、Overfull box はない。
既存の小型数式フォントのサイズ代替警告は残る。
新節をテキスト抽出で確認し、非平方 a=0 の表・証明を含む主原稿101ページと
評価追補2ページを画像でも点検した。
最終ログは [主原稿](MANUSCRIPT_SYNC_20261103/main-pass3.log) と
[評価追補](MANUSCRIPT_SYNC_20261103/status-pass3.log) に保存している。

原稿索引の v61 の3主張IDは `research/CLAIMS.json` に存在し、
v60/v61 の全参照ノートのファイル存在を確認した。数学証明コードをこの原稿同期で再実行したわけではない。
