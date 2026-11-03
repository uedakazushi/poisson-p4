# CA/CB の局所正規化・交叉に関する独立原稿

2026-11-03。新しい日本語 LuaLaTeX 原稿
`papers/ca_cb_local_normalization_ja.tex` と PDF を作成した。
最終 PDF は実測 **17 頁**である。研究登録と main への保存は統合担当が行う。

原稿は共通階数二点での個別の被約成分 CA、CB の全形式完成型、conductor と二重被覆の
rigidity、新しい Cohen 底からの pushout、完全な 99 関係、seminormalization の完全な
81 関係を扱う。形式的な特異・非正規・非 CM・非 seminormal locus と正規化 fiber も
証明する。さらに実際の二つの多項式自己旗族 Pa、Qb を使い、非空 open 上の actual
scheme intersection D における二つの異なる全完成型と、射影 19 次元の global embedded
associated component の閉包を扱う。Pa の nil module の generic 長は 1、Qb は 4 であり、
Qb の embedded primary 商長 6 と区別している。

これらの embedded component は全 Poisson scheme の新しい被約既約成分ではない。
全分類の下界 58 は変わらない。共通基点の相対的な全 CA/CB 交叉、global exhaustion、
収束する正則同値、文献上の優先性は主張しない。個別の形式座標変換を相対配置に流用しない。
全非零 parameter での旗数と、unit 条件を要する全 D 完成型の非空 open は分けて記した。

有限入力は元の全 150 係数 source、全 proper flag fiber、実際の odd 関数の jet、
全関係 ideal、独立監査の受理済み report に結び付けた。とくに Qb の独立三次監査の
`independent-source-curves.json.gz` を直接保持する。新しい checkout での cache 復元 helper
`code/restore_CB_cubic_worker_inputs_20261103_v63.py` と、保存証明書の監査および明示的な
三次 worker 再計算の違いを付録に記した。

Cloud の Macaulay2 1.24.11 による default GB の全 99 関係 kernel 等号は、保存済みの
独立証明書を参照した。この原稿作成で重い CAS 計算を再実行していない。関係の最小性も
主張しない。LuaLaTeX 最終 pass に未解決の参照・引用、再実行要求、overfull box はなく、
1 頁と 16 頁を画像として確認した。日本語 italic 情報と小さな数式 font の置換は build log
に明示して残した。PDF は A4、form なし、JavaScript なしである。

再現可能な保存入口は
`code/ca_cb_local_paper_manifest_20261103_v63.py`、結果は
`reports/CA_CB_LOCAL_NORMALIZATION_PAPER_20261103_v63/` にある。
そこに source/PDF SHA256、実測ページ数、最終 build log、入力ハッシュ、独立証明書の
検査結果、登録依頼を保存する。従来の 118・16・91・14 頁の各原稿と PDF、および
歴史的な `papers/SOURCE_INDEX.json` の合計 9 ファイルは作成前の SHA256 と一致する。
外部投稿・公開サービスへの送信は行っていない。
