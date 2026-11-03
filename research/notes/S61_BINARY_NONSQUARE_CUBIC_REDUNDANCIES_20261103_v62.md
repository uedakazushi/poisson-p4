# 元135式の全三次matrixから三つの重複rowを除く

2026-11-03 / v62。既存の63×91全三次matrixについて、core開 a*d*F*B*C*J≠0 の全点で原row48,14,55を順に削除できる三つの厳密なQQ多項式恒等式を得た。番号は0始まりである。残る60×91のrow moduleは元と同じになる。これは各点の階数57下界や四次の全被覆ではない。

ここで F=500a+90e²+3e−24、B=(50a+18e²−9e−2)d+60e−40、C=100a+18e²−9e−2、J=5ad−6e−4。Jは他のframe因子を指す記号ではない。

入力は [現在の135式から再構成した全三次matrix](S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62.md)。そのfull入力のSHA256は `dae59904925260605628d925ae0565de084c9af0ae1bccf1774f1af067830ad6` である。今回の検証器は凍結archiveからこの入力を読み、hashを照合する。

各rowをa,F,B,C,J,dの共通因子と非零QQ定数だけで正規化し、元rowが正規化rowと記録した単元倍率の積に等しいことを全成分で検査する。正規化後のrowをR_iとする。三関係の係数次数は順に2,3,6であり、削除するrowの係数はそれぞれ1、12a/5、1である。後の関係は前に削除したrowを含まない。各関係を全91列に掛けた結果は、QQ[a,e,d]の多項式として厳密に零となる。

従って各関係から指定rowを残りのrowの組合せとして解ける。逆に残りは元matrixのrowなので、row moduleは両方向に等しい。係数次数6はcore係数の次数を指し、全matrixの法線重み三と混同しない。

関係の発見には有限係数次数のMacaulay行列とsparse有限体演算を使った。有限体rankの下界、rational reconstructionした候補、異なるfree coordinateでのidentity、元matrix全列でのQQ恒等式を合わせるrank sandwichにより、固定した次数以下のkernelを検証した。全syzygy moduleの生成は主張しない。公開の [QQ検証器](../../code/s61_binary_nonsquare_cubic_redundancies_20261103_v62.py) は、その探索をやり直さず、保存した有限係数データを元入力から独立に再検査する。

一般階数57と全点上界57は [既存の幾何学的証明](S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102.md) を使う。今回の三row削除後にも一般依存関係はさらに三つ残る。60行全てが独立であるとも、57行への全開での縮約を完了したとも主張しない。

再現:

    source /workspace/cloud-setup/activate.sh
    .venv/bin/python code/s61_binary_nonsquare_cubic_redundancies_20261103_v62.py \
      --output /tmp/nonsquare-cubic-redundancies-new

凍結 [report](../../reports/S61_BINARY_NONSQUARE_CUBIC_REDUNDANCIES_20261103_v62/report.json) に全入力・生成hash、単元正規化の記録、全QQ恒等式のnative log、残るmatrixを保存した。V4の探索・負の有限次数結果・P=adの別座標は [引継ぎpacket](../handoff/20261103_v62/v4_sparse_kernel/README.md) に分けて保存している。
