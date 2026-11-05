# col53 の完全 H 二面：符号付き値の独立検証

独立検証は PASS した。報告 SHA256 は
`d1708919e26ff8d484855c987d432cc4e5245cab8ae58164c5862bb44a348c26`。
元の SOURCE/code peer `7c6e2587…` に続き、実値 primary `b8d4309e…` を固定した。
元の1,203 SSA nodeから H=−3 と H=21 の二面だけを再計算し、
保存された実値との全係数比較に成功した。

H=−3 の正 SOURCE 支持は766項であり、符号付き実支持は411項となった。
有理数 content は `1/503316480000000`、全係数 hash は
`4813aeeda3b145c587ed8c551b6a5d51d188f2e17701f429facb175a4ed73688`。
H=21 の正 SOURCE 支持は214項だが、符号付き面は厳密に零である。
九 byte の零 cell 表現まで別 decoder で確認した。

独立実装は author の producer・native・codec を import していない。
既存の独立 vertex reviewer の literal selector と primitive binary decoder だけを再利用し、
入力項・積の二重 loop・和の子の順序を逆順にして整数 map を再構築した。
全五座標を保持し、source q と既知 QQ/monomial 変換を固定して計算した。
キャッシュ cell の九 byte header と、UF10/C10 の係数 `353894400000`、
shift `(1,0,0,0,1)` を厳密に照合した。

面の枝刈りは source q の厳密な零および source H の一致だけによる。
先に零と分かった十個の lex 頂点係数も、保存済みの実面の相殺も利用していない。
独立 root map が完成してから実値 cell を読み、411個すべてについて係数・五座標・
primitive gcd・content・正規順序・全値 hash を照合した。

具体的に指数 `(-3,0,6,-1,10)` の係数は
\[
\frac{241815}{134217728}\ne0
\]
である。したがって固定 row6/col53 の分母を払った Laurent 分子は非零である。
これは全1,203 nodeの最終式の一面の非零であり、411項の面を全式と同一視しない。
ほかの H 面、全成分 map、階数・新しい pivot・全 inverse・尽くしについての主張は加えない。

独立 reviewer の metadata charge は102,267、事前の全所有上界は130,358、
実所有最大は128,188で hard150,000内だった。五座標項ごとに7 slotsを課した。
leaf decode16,152と同時存在 mapを含み、全 unique RAW入力は14,738,073 bytesである。
before-operation 整数上界は実測174 bit、保存整数は173 bit。
committed640・temporary1921・hard2048を変更していない。
55秒 signal、2 GBと途中 prefix 保存を備え、実行は1.209秒で終了した。
新しい parameter point・prime・元係数候補は追加していない。

再現コードは `code/s61_col53_H_faces_values_review_20261105_v64.py`。
凍結 report と primary は変更せず、再実行が必要なら別作業コピー・別出力先を使う。
