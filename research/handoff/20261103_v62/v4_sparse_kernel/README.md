# V4 sparse kernel / P=ad chart の保存状態

2026-11-03。親agentから service capacity error の報告を受けた後、最初に現状を固定した記録である。モデル変更・成功済み計算の再実行はしていない。V3の保存済み公開物は変更していない。これは研究の停止指示ではなく、引継ぎ先はこの状態から続行する。

最後の成功は P=a*d の係数座標で、三行を削除した60×91 matrixに対する係数次数≤6の有限Macaulay行列が rank5040、kernel0となったこと。新しい行削除はない。native exit0を確認した。保存時点で自分の Singular / kernel solver の生存processは **0件**、live PIDは無い。cgroup memoryは約5.84GB、oom/oom_killは0。親の他の研究processを停止していない。

正の追加成果は、元core座標の61×91 matrix（原row48,14を削除済み）で、係数次数≤6の5124未知数のkernelが正確に1次元と証明し、全QQ多項式恒等式を再検査したこと。関係の原row55の係数は1なので、この行も指定開a*d*F*B*C*J≠0の全点で削除できる。現在は60×91、一般階数57が目標で、さらに三つの一般依存関係が残る。全点での階数57下界、四次の全被覆、全非平方層のambient除去は証明していない。

## 厳密な計算の意味

`SPARSE_FINITE_MATRIX` は、core係数次数を固定した左syzygyの有限係数行列を表す。係数体はQQ、native discoveryは素数2^61−1を使う。sparse echelonから rank の下界を得て、rational reconstructionしたkernel候補の**元matrix全列でのQQ多項式恒等式**と、異なるfree coordinate上のidentityを検査する。これが成功したときだけ rank sandwich によりその次数以下の全kernelを証明する。全syzygy moduleの生成を主張しない。

`early` 版は、全入力equationを処理し終える前でも、必要な個数の独立QQ kernelが全式を満たせば同じrank sandwichで終了できる。dense RREFの一方向の候補検査や、非斉次modStdのraw unit出力を採用する方法ではない。小fixtureと、既存のexact QQ次数2kernelで独立に照合した。

係数次数の負の結果は、**それぞれの固定matrixと係数座標**に限る。元座標では、二行削除後の次数≤4,5、三行削除後の次数≤7で追加kernel0。P座標では三行削除・既知単元正規化後の次数≤3,4,5,6でkernel0。従って、残る全syzygyが無いという結論ではない。

## P座標

第三core slotを P=a*old_d とする。aは既に単元なので可逆な座標表示である。各単項式は a^i e^j old_d^k → a^(i−k)e^j P^kとなり、各rowの分母をaの冪だけで払う。全変換は `product-core-chart.py`、単元倍率とhashは `P-chart.json` にある。native形式の第三slotの文字名はdだが、ここではPを表す。

単元は a,F,B0=a*old_B,C,J=5P−6e−4,P。必要なpolynomialは `P-chart-units`。a,F,B,C,Jの未確認の別因子を逆にしていない。

## 再開点

原matrixと原正規化matrixは保存済みV3にある。V4の大きい導出matrixはmanifestへhashを記録し、bundleには重複収録しない。`make-early.py`で原row14,48,55を削除したmatrixを作れ、`product-core-chart.py`とnativeの単元正規化からP matrixを再生成できる。開始commitはroot側で記録済み。この保存時のmainは0713e279a55bcf81590242e3b908509a3cefa499。

次は、残る三つの一般依存関係をP座標の次数7以上、または直接のJacobi恒等式から得る。既存の成功済み低次数計算をやり直す必要はない。native sparse solverの可変部分は `sparse-early-certified-kernel.cpp`、CLIは `binary plain-matrix core-degree output [known-unit-file]`。compileは既存activate.shを読み、`c++ -O2 -std=c++17 source.cpp $(pkg-config --cflags --libs flint)`。

`code/data/s61_binary_nonsquare_cubic_redundancies_20261103_v62.json` は検査済み三関係の有限データである。新たに書いた `code/s61_binary_nonsquare_cubic_redundancies_20261103_v62.cpp` は、それらを独立に再検査する verifier の**未実行draft**であり、保存時点ではまだcompile/replayしていない。これを実行済みの証明書と読み替えない。公開用Python producer/noteはまだ作っていない。reviewer connectivityがV4のbounded reviewを担当する。
