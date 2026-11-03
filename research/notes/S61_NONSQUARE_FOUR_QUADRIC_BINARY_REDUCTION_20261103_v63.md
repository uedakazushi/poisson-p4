# 四二次式の四線形変数構造と pair gate への縮約（2026-11-03）

実際の804係数項から再構成した式に対し、`W=x5`, `c=6075`、既存の `A,B,Q2prime` を用いる。QQ 恒等式

`Q1=AB/c+W*C1`, `Q2prime=40*x0*B/c+W*C2`, `D=40*x0*C1-A*C2`

が成り立ち、C1,C2 は線形式である。従って Q1,Q2prime,D/c は線形行列 `[[A/c,40*x0/c,W],[C1,C2,-B]]` の符号を変えた2-minor である。

6個の線形 entries の rank は、想定された6ではなく、**全 a*P*J!=0 上で正確に4** である。実際に `(A,x0,W,B,S,U)` を座標として、逆変換の行列式が `1/(200*a*P^3*J)` となることを検証した。この変換後 C1,C2 は S,U を全く含まず、A,x0,W,B のみの線形式になる。すべての係数と逆変換は certificate.json に保存した。未知因子を逆にしていない。

その結果、B=W=0 上の厳密な式は

`D = -8*(3Pe+2P+12a)/(1215P²)*A² -256a²J/(243P)*A*x0 +1280a³J²/243*x0²`。

A と x0 はこの制限空間上で独立である。例えば B=0 で `A=(10aPJ/3)*x4`、B には x1,x3 の非零係数があり、`3Sc-7Vc=20P` がそれを保証する。従って D のこの制限は全域で零でなく、表示した mixed coefficient と x0² coefficient は既知 unit である。

**pair gate が成立する点では、補助五二次式の rank5 も自動的に成立する。** ここで pair gate の制限空間は `A=x0=W=0` であり、**B は自由** である。五二次式の任意の定数係数関係をまずこの空間に制限する。Q1,Q2prime,D は零であり、pair gate によって Q0,Q3 の係数が零になる。次に残った三式の関係を B=W=0 に制限すると上の非零 D によって D の係数が零になる。最後に W=0 では関係が `B*(c1*A+40*c2*x0)=0` となり、B≠0 と A,x0 の独立性から全係数が零となる。

既存の gcd gate と線形 syzygy の縮約補題に組み合わせると、gcd gate と pair gate が成立する点では raw 24×56 の rank24 が従う。gcd gate の例外全体は別の finite Nakayama packet で actual raw QQ surjectivity が既に証明された。従って残る十分条件は pair gate の全域証明だけであり、別の五二次式 rank gate を追加で解く必要はない。これは pair gate 自体の証明ではなく、全 QQ product41/full57/ambient removal を主張しない。

再現: archive 展開後 `python code/s61_nonsquare_four_quadric_binary_reduction_20261103_v63.py`。約2.6秒で全804項、逆変換、三つのminor、全線形係数、D の unit係数を厳密QQ演算で検証する。GB、有限体、数値サンプルを用いない。
