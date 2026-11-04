# native128 X1：独立C36・実source・off-intersection監査

2026-11-04。作者workerやhelperを呼ばず、先行の独立C4 tensor engineと疎dual演算を再利用し、
別の辞書基底C36積を実装する。`epsilon_i epsilon_j=0`、`tau_i²=nu_i,1`を文字どおり使う。
tauのcarryがbase一次係数を作る項を保持し、両tau二乗を最初から零とはしない。
maximal ideal五乗零の有限geometric inverseと、定数有理逆からの四段nil摂動で行列逆を作る。

既存105解のArtin20値から同じCA-F9八座標の一次係数を抽出し、全17926 source DAGを疎dual9で
再評価する。物理150 trace pairとnuの全一次値を照合する。source48/105を解き直さず、
新しい物理source jetを生成しない。

保存nil中心の同じ十五constraintに対するlength4での固定225 mixed係数を独立tensor pushで照合する。
J30のblockと保存逆の両積、240 RHS、240補正係数を確認する。保存補正をC36へ代入して全30成分の
constraint零を直接確認するので、未知のgauge解を仮定しない。J30はunitであり、このpointed解は一意である。

全finite flag/source式、physicalF12、R17を再構築する。C36ではF12は実際に非零である。
`F17=R17 F12`を確認し、最初五pivotには完全なordered polynomial telescopesを使う。
原点でだけ正しいJacobian shortcutへ置換しない。全R12、T、N、G8を照合し、

\[
 dH12=R12F12,\qquad TH12=(G8,0),\qquad NF12=d(G8,0)
\]

を全36成分で検査する。dを反転せず、`dG8=0`からG8零を推論しない。
全constant N/G8/Sは独立受理済みC4 origin peerに一致する。
最後に正しい転置方向の`Res4(N^t S)`を独立乗法で再構築し、全128列・6144個の定数とX1を比較する。

元の係数関数への正式な第一微分であり、I²や348環への先行reductionは行わない。
X1のrank、T最小生成数・長さ、K torsion、元深さはこのpeerでは計算しない。
専用checker `code/ambient_B18_native128_first_review_20261104_v63.py` と同名REVIEWに保存する。
