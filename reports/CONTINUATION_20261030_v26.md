# v26-s61 研究checkpoint

2026-10-30（Asia/Tokyo）。開始commitは `59ddf2cdd92c338fdc4006f2d48dae9c820189d6`。
統合先はmain。成分数の下界57は変更せず、完全分類は未達成。

## 数学的な進展

1. T01のprimitive一次部分への射影の優越性。6個のtoricパラメータとSp4の10方向の
   微分に、非零16×16小行列式 `258241/10504375` がある。
   一般一次部分は既知toric構造へ延長できる。T01成分上では一般に有限、次数は未決。
2. 全S61の座標adapterと、toric6変数・almost-abelian2変数の一般有理写像を構成。
   両族で150復元係数、134式、projective Jacobiを厳密に検証。
3. 有限perspective公式から自然なprimitive gaugeを構成。
   非零係数階数1は5次元の単一Sp4軌道であり、その全解は固定fibreへ帰着。
   固定fibreは45変数・134二次式。既知点の接空間次元38。成分分解は未完。

## 実行と再現

- `bash scripts/setup.sh`：成功、Python3.12.14、固定依存版一致。
- `.venv/bin/python scripts/verify.py --suite all --jobs 4`：integrity、1298項目、10群すべて成功。
- S61原コードの再生成：5 JSONはruntime以外完全一致、Singularはバイト一致。
  復元後Schoutenで154式と55消去行を独立比較。`S61_REBUILD_20261030.json`。
- 新座標adapter：24検証群成功。符号、Euler gauge、有限群作用、分母付き恒等式、負例を含む。
- 一般族：`S61_FAMILIES_20261030/verification.json`。有限体・数値標本に置換していない。
- primitive gauge：17項目。rank1表現・fibre：有理数体上の行列式と全式を保存。
- 優越性証明書は別担当が再実行し、統合時にも保存行列から行列式を独立に再計算。
  `S61_INDEPENDENT_REVIEW_20261030.json`。
- JSONと依存ID・参照先、主要コードの入力SHA、Python構文、差分空白検査を確認。

一般族の全式代入では、汎用の有理式展開が遅かったため、共通分母を払う疎多項式比較へ
変更した。途中停止した旧方式を成功とは記録せず、新方式で全134式が完了した。
階数1fibreでは134式の定数線形包に非零一次式がなく、その方法だけの変数削減はできない。
この障害はイデアル全体の一次式の不存在や計算不可能性を意味しない。

## 保存と環境上の未完事項

原資料は変更していない。新原稿は日本語LuaLaTeX、
`papers/s61_first_jet_reduction_ja.tex`。ltjsarticle、luaotfload、HaranoAjiが不足しているため
今回のPDFは未生成。新規性の一次文献監査は未実施で、内部草稿として保持する。

既存bundle生成コードに旧BUNDLE_MANIFEST.jsonの重複格納があったため、今回のmanifest
だけを生成するよう修正した。bundleはCRC・重複エントリ・全ファイルハッシュを確認する。
mainへ通常のcommit/pushで保存し、継続ZIPをユーザー向け出力へ作る。
独立したLibraryアップロード機能はこの実行環境になく、Libraryへの保存は未実施。
Pagesのファイル添付は10MiB上限で、この原資料込みbundleの代替にはならない。
成功した保存先とcommitは最終応答およびGit履歴で確認できる。

## 次の実行可能な課題

`research/NEXT_TASKS.md` を入口にする。T01の優越性を前提に、一般一次部分上の
全延長、他成分の一次部分の像、または階数1固定fibreの分解へ進む。
まず追加成分CAまたはCBの一般族を現行adapterへ写し、primitive一次部分への
写像の微分階数を求める。rank1fibreを選ぶ場合は、保存45変数134式から安定化群や
局所化を使う削減を明示する。大規模primary decompositionを理由なく再実行しない。
全51成分の個別S61カタログ、全体の根基分解、特殊軌道・交叉・冪零構造は残る。
