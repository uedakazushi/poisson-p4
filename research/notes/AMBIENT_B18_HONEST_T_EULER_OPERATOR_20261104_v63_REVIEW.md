# 実 T の116生成936関係：独立 source・Euler・完全性監査

2026-11-04。作者の producer は import・再実行しない。既受理の実
`P:640→9`、`UP=QB`、完全 `Syz(Q):128→48`、比較単位 section、
実 Cbar の三 sector Koszul 解像度を hash で固定する。

独立 checker は最初640列の `[B48;P8]` と後128列の
`[Res4(N^t) Res4(S12×32);0]` を各原 source 参照まで再構築する。
比較行列の符号は `V=[Q48,-U8]` なので、前者は `QB-UP=0` である。
`U(1)=1` の単位 source と定数 target を消す際、`UP=QB` の定数行が
P0 を正確に与え、残る56 source 座標や768列の符号は変わらない。
旧 regular r/g は元の正則関数として制限し、特異商で割り直さない。
P の640関係の完全性は、別の受理済み honest C0 peer を明示的に追加して固定する。

選択24列 M の原点逆行列は両積を独自 Fraction 演算で検査する。
その formal inverse の行列入力と全 n の参照先も確認する。
source 同型は `(a,b)→(Ma+Nb,b)`、逆は
`(c,b)→(M^-1(c-Nb),b)` であり、`d1 m=V` から V は `[d1,0_32]` になる。
したがって核は `im(d2)⊕Rbar^32` である。

三 regular sequence は、実 CA-F9 の同じ r と z 七座標で
`(z,r³)`、`(z+c r²,r³)`、`(z,r²)`。中央だけ `u=z+c r²` に変える。
全84 d2 と168 d3 を exterior basis と符号から独立再構築し、作者が保存した
全70017 coefficient slot、74176 typed node、全768列の Euler readback、
補集合32座標の付加を照合する。単なるサイズ一致を source の同一性としない。

Euler 公式の全次数証明は次の通り。u-次数 d の v_j と次数 d+1 の v8 に対し、
cycle `Σu_i v_i+r^n v8=0` を微分すると

\[
 \sum_i u_i(\partial_i v_j-\partial_j v_i)
 =(d+1)v_j+r^n\partial_jv_8.
\]

微分後の次数に `(E+2)^-1` を適用するので左の wedge 係数の分母は d+1。
`w_j8=(E+1)^-1∂jv8` も分母 d+1 であり、Koszul differential の
`-r^n w_j8` が追加項を消す。最後の成分は各正の u-次数 k ごとに
`Σu_i∂i v8/k=v8`。次数零部分は cycle を u=0 に制限して
`r^n v8(0,r)=0`、正則 domain `QQ[[r]]` により零となる。
よって全形式級数に `d2 w=v` が成立する。

ここで分母は正の整数だけで、r や非単位 source 関数を反転しない。
微分は次数を高々一つ下げ、Euler 逆は係数を次数別に再尺度化するので連続である。
中央の triangular shift とその逆も連続である。従って各実係数関数の全級数に作用する。
有限多項式 fixture は符号確認だけに使い、全次数の証明の代用にはしない。
checker は別の係数 scatter 法で、作者と異なる高い u-次数および r^4 sector も確認する。

contraction 自体は Rbar 線形ではない。768個の固定列それぞれに preimage を選び、
その後にその列から自由加群写像を Rbar 線形に延長する。
`Rbar^84⊕Rbar^32→ker V` の核は `im d3⊕0`。
各 X 列の選択 preimage とこの核を合わせると `im X` の全逆像になるため、
168+768=936関係の cokernel が実 T と正確に同一である。
未知の Syz(P) を落とす操作や、有限 fixture からの外挿はない。

操作は**元の適合形式関数を I² で割る前**に行う。Euler 逆から得る係数は
形式級数であり、代数関数であるという主張はしない。
X(0)、X の一次、極小化、348環の係数、T 長、K torsion、元 B18 深さは未計算。

再現は新しい保存先を指定して
`PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_honest_T_Euler_operator_review_20261104_v63.py --output work/new-honest-T-review`。
作者の凍結ファイルは変更しない。
