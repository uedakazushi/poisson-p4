# 実比較X(0)：独立C4 source・adapter・固定rank監査

2026-11-04。作者のproducer/helperを呼ばず、元の有限source式を
`C4=QQ[tauA,tauB]/(tauA²,tauB²)`で独立評価する。新implicit solve、source jet、
点、prime、GB、X1、348係数行列は計算しない。

四成分演算は閉じた積公式とunit逆公式を使う。行列逆は、原点有理逆Cとnil摂動Eに対し
`(I-E+E²)C`で作り、E³零と左右逆を確認する。作者のC4 Gaussianは使わない。
HA/HBの各source first direction、CAのactual group15 gaugeを同じraw150基底で組み、
外積二乗Gと対称二乗G逆のliteral tensor pushを行う。作者のLie helperを用いない。
全14 source値、相対五flag、全17式、実retained κ9を照合する。

同じfull150作用からR17を再構築する。最初の五unit root、実Jacobian Schur、row unitから
R12を作る。次の五unit rootとordered prefix synthetic divisionでTを作り、実CI graph cofactors、
Q0/Q1/sr²の正規化、全12 source rowと照合する。N=T R12と保存逆の両積を確認する。
いずれの段階でもdや消えるsource関数を反転しない。

特に実八生成元は

\[
G_8=(1/2,25/8,0,0,0,5,-125/2,0)\,\tau_A\tau_B
\]

であり零ではない。`d G8=0`をG8零と取り違えない。このベクトルのKoszul mapは
QQ係数rank7、`tauA*tauB`乗算はrank1なので、28 native列のRes4 rankは7。
四つのzero-source自由列のRes4 rank16が加わり23となる。
unit N転置がこのrankを保つ。転置方向と元のnative S符号を保持し、全48×128成分を比較する。

最初640列のBとhonest P8は、先行独立監査済みのregular origin cacheから全て評価する。
旧r/gを零locusで割り直さず、元の正則関数node IDを使う。全768列のm0 X0零を確認し、
作者が保存した固定消去の引算・正規化・元source線形結合を全768列で再検査する。
異なるpivot行のunit triangular性が下界、残る列の明示的零化が上界を与えるため、
新しいrank searchをせずrank27を確定する。全X0は56×768、非零係数97個である。

これだけから得るTの最小生成数下界は32−27=5。X1に由来するEuler wedgeの定数が未計算なので、
Tの最小生成数・長さ、K torsion、元B18深さを決めない。
専用checker `code/ambient_B18_comparison_X_origin_review_20261104_v63.py` と同名REVIEWに保存する。
