# col53 十頂点係数の独立 SOURCE・VALUE 検証

独立検証は PASS した。報告 SHA256 は
`d7e8134bbed2dc1483bef0933d9c892223222583619f72583d564ea2dbb66831`。
固定 row6/col53 の cleared numerator について、指定された十個の普通の
Laurent 係数は全て厳密に零である。成分全体の零・非零、pivot、rank、
極値面全体、補集合の尽くし上げは結論していない。

原 SOURCE gate は `6ad0a98af6dd6539bf2d9cdf16ec8e431ee70df66dbafef540c6aac32c73a1b8`、
原 VALUE report は `ef9cb1c8175ce59df1cd95274bcc92df719036a61b4d7d41d080237e4b25ea20`。
producer・native codec は import せず、元八 phase、45項と48項の計101 receiptから
SSA を再構成した。全1,203 node の演算・親参照・有理数 scalar・単項式 shift・
元 source binding が一致する。raw53 の各再定義は独立した SSA ID として扱った。
fused union の literal_dot が元 A,B の積であることも確認した。

入力91本は個別に元の literal JSON または binary payload へ戻った。
整数 varint・符号長・primitive 有理数 content の codec を別実装し、
offset 付き payload と Q5CELL1\n＋flag の九 byte header を持つ cell を区別した。
col56 の UF10/C10_Ni56 は scalar353894400000 と shift(1,0,0,0,1)、
他の既存 cell は scalar1・shift0 であることを明示的に確認した。
JSON binding の元 literal と QQ/monomial 変換も保持した。

各 source 支持に対し、五 cyclic lex 順序とその逆順を固定する。
これらは Z⁵ の平行移動不変な全順序であり、有限 Laurent 支持上で最小元を持つ。
子の全単項式はその正支持最小元以上にある。積で二つの非負な差の和が零になるのは
両方零の時だけなので、親の正支持極値の係数は子の極値係数の積となる。
和の正支持極値に寄与するのは、同じ極値を持つ子だけである。

ここで「正支持」は係数の絶対値上界に使った source 支持である。
子の極値係数が相殺して零になっても、その元の正支持極値を保持する。
actual の次の頂点に置き換えないので、この帰納計算は相殺時にも正確である。
これは固定された十個の係数 functional の計算を証明するもので、
零になった係数から子や根の多項式全体の零を証明するものではない。

全12,030個の node/order 頂点座標を独立に再構成し、最終 SOURCE snapshot の
根の全36,620キーを stream して十座標を照合した。
その後、事前 SOURCE/budget receipt を保存してから、順序一つずつ独立に
有理数の係数を伝播した。1,203×10の係数配列を同時には保持していない。
得られた十係数は原 VALUE の全十記録と一致した。

入力・node の正値 profile 上界402 bit、既存 committed640 bit、
全一時演算1921 bit／hard2048を区別した。各積・和の**演算前**に
分子分母の bit 上界を検査した。実際の最大は237 bitであった。
独立な所有数計測は、JSON 後の list に入る全指数ラベルも数え、
tuple/list の事前計数差に依存しない。正値の全所有数上界118,477、
RAW入力と予備は32 MiB未満で、60秒・2 GiB・150,000 slots以内だった。
実行時間は0.894秒である。

初回 CELL header reader の失敗コード・ログ・receipt も不変のまま保持し、
今回の V2 の source と変換を別に照合した。旧 phase・原報告・原 VALUE は変更していない。
全167×14、87追加 RHS、66制約、one P290、F9/F10と全補集合を保持する。
新しい parameter point、prime、全成分展開、rank search、極値面 VALUE は実行していない。

再現 checker:

```sh
.venv/bin/python code/s61_col53_vertex_source_values_review_20261105_v64.py
```

凍結出力先への再実行は拒否する。再検証は既存成果物を保存した別作業コピーで行う。
