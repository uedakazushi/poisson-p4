# 次の課題 — 全被約成分の尽くし

2026-11-06 / v64-classification-first。classification/PLAN.jsonが現行キュー。

```sh
.venv/bin/python scripts/classification.py check
.venv/bin/python scripts/classification.py next
```

C0とE0の有限準備は完了。主課題Gは未解決。
まずclassification/GENERIC_NORMAL_CHART.mdと独立reviewを読む。
既存splitter A=E5×Bを用い、tauを超越のまま保つell=Σtau^j n_jにより、
元134式+1-u ellの一つの単位証明書で全未知残余を検出できる。
同じK5源の259法線式も全schemeを保つが、係数体solveの展開を先に必須化しない。

保存済み中心GF(32003)(tau)のSymPy F5bは55 CPU秒で未完了。
同じbackendの点・次数だけを変える探索を繰り返さず、native exact backendまたは
構造的syzygyを選ぶ。成功時にもQQ/一般底への正当な移送と底分母の独立照合が必要。
固定normal0の診断も305wall秒で未完了であり、直接SymPyGBの反復は停止。
全K16上の固定eta一つの135式判定はSYMMETRY_FIXED_NORMAL_CHART.mdにある。
42軌道spanは独立証明済みだが、K5切片への非不変形式の輸送は未証明。群変数を無断で落とさない。
現在は標準network権限によってSingular4.4.1をwork配下に用意できた。
Native stdの300CPU/315wall/2GBジョブも未完了だった。次は既存G6帰結の厳密前処理とslimgbの別SOURCEを点検する。
公開パッケージhashと起動環境はreports/CLASSIFICATION_NATIVE_ENV_20261106.json。

Eではclassification/EXCEPTIONAL_BASE.mdの必要条件を使う。
一般五点の底開D(f_G)なら、f_Gは既知49特殊閉像すべてで零となりG28を因子に持つ。
この条件は十分ではない。底の判別式だけで特殊部分を尽くさない。
beta=0の全134行は45未知数・70二次式+64零行として保存・独立確認済み。
旧Q45の22族から既知52族閉包への全包含はBETA_ZERO_COVERAGE.mdとE_ZERO_REVIEW.mdで閉じた。
次は他の特殊底で、同じ被覆論理へ接続できる最大次元の層を選ぶ。

Aではclassification/KNOWN_COMPONENT_AUDIT.mdの弱いrank2接続を使う。
八rank2-locus族すべての成分性を主分類の前提にせず、六全空間成分と
closure(rank4)による支持被覆を確認する。安定接層証明・追加11族とNEW25の
成分性・相違・一次引用の残る荷重を主原稿へ接続する。

主原稿はpapers/poisson_p4_classification_ja.tex。
PDFは現TeXから8頁で再生成・確認した。G/E/Aが閉じるまで分類主定理を未証明と表示する。
B18 H0/T2/nil指数、H絶対既約性、別論文の仕上げは休止のまま。
checkpointは完了宣言ではない。利用可能な証明・厳密計算で実研究を続ける。
