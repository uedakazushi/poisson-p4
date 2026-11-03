# 全35変数の真の三次初期部分：独立監査

2026年11月3日。全350式の QQ 三次 DAG 抽出と全150係数からの physical350 Jacobi を再実行し、
保存 archive と一致した。さらに全重み block の GF rank と、分母を払った整数行列の FLINT QQ rank を
再実行し、1894→1911、新17 classes、通常 Hilbert function (1,35,573,5859) が一致した。
抽出14.95秒・QQ native32.04秒であり、大きな GB は再実行していない。

全350行の保存された定数 quadratic elimination と cubic carry を Fraction で独立に再生した。
57 pivot と293定数 syzygy、135非零 cubic initial が正確に得られる。
真の I3 がこれらと全1995本の linear quadratic multiples の和で尽くされる証明も受理する。
一般の ΣhiFi の三次部分は hi の定数・一次部分だけに依存し、定数の二次消失条件は保存した
293 syzygies が尽くすためである。単に135本の裸の cubic terms を数えていない。

旧13座標平面への全制限を独立に計算した。非零は source234 の3/16、source349 の−1/16倍の

    z8²(12z14+z17+6z3)

だけである。従って次数≤3で生成した ideal の **この平面への制限** は厳密にこの主 ideal。
二つの12平面の generic multiplicity は2と1である。
また quadratic の square monomial が禁止しない15変数について全32768部分集合を独立に調べ、
最大 coordinate plane13が旧平面だけであることを確認した。z8=0 の coordinate12平面が残るので、
J≤3 に入る coordinate plane の最大次元は正確に12である。

追加の17本の明示 cubic representatives の frozen native selection も再実行して一致した。
QQ rank差17と同じGF rank、unit pivots によりこれらは QQ と Z_(32003) で quotient の基底を与える。
次数2・3の有限線形代数の image・quotient が free になる根拠は Smith rank と residue rank の一致であり、
全局所環の flatness や higher initial ideal の base change をそこから推論しない。

受理するのは真の I2,I3 と J≤3 の明示57+17生成、旧平面の制限までである。
全 tangent cone の Krull 次元≤12、実12次元枝、全 completed ideal の74生成、higher obstruction の消失は
いずれも結論していない。code/q_ordinary_cubic_review_20261103_v63.py と同名 REVIEW report が証拠を保存する。
