# K7W 固定9次minorとQQ gcd frontier：独立監査

2026-11-04。新しい source や target、点、prime、minor、NF、rank、因子を探索しない。
既受理の K7W 一行 peer、pure90 の元M15×98参照、既存五minorのH15 gateを固定する。
作者の producer や helper は import・再実行しない。

固定行列は旧case9の行 `[9,2,3,5,7,6,1,4]` と、literal normalized C68 と
physical v3=joint W の積である原319行226を第九行にする。列は
local98の `[0,1,2,3,4,5,6,9,28]`、original122の
`[0,1,2,4,5,6,7,10,33]`。独立 reader は九つの元行ファイルから全81値を取り、
既知五曲線unitの各row LCD、1057個のRF入力項、元flat inputの全byteを照合する。
source normalizationと全210/NF98は先行K7W peerを再利用する。

独立native checker は81個の分母が各既知unit LCDを割ることと clearing を正確に検査する。
作者のBareissとは異なり、固定列集合ごとのLaplace展開を隣接二層だけ保持して計算する。
この計算は行列のpivot・逆やrankの探索を含まない。
保存されたN202と全係数で一致し、`N=H Q+R`、次数15/187/14も文字どおり確認する。
元QQ入力とcleared係数、H,N,Q,Rの分母は32003-unitである。
NやRのleading coefficientのunit性は仮定しない。

QQ gcdを再計算する代わりに、保存済みの二つの非零scalar contentとprimitive多項式を使い

\[
 H=c_H H_p,\quad R=c_R R_p,\quad
 g=S H_p+T R_p,\quad H=g Q_H,\quad R=g Q_R
\]

の全QQ係数を独立に検査する。monic

\[
g=z^2+\frac{200}{13}z+\frac{22800}{13}
\]

がこの両方向の恒等式を満たすため、厳密にQQ gcdである。
mod32003の次数だけからQQ gcdを推定していない。

旧五minorについてはそのsource-bound determinant gateを受理済み前提とし、保存された
五つのprimitive/H除算と、p-integralかつleading coefficientを保つquotientの
既存prime cofactor Bézoutを独立に再検査する。これにより同じH15の使用を確認する。
既知unitをstripした旧五minor ideal `(H)` に新しいminorNを加えると、許容曲線unit ringで
`(H,N)=(H,R)=(g)`。この**選択した六minor**の共通frontierが二次有限QQ代数と定まる。

従って既存の全P100 pointwise149を保ち、`D(g)`では150下界を与える。
gの零点における実行列rank不良、全曲線150、full239、global9target inverseは主張しない。
作者のcompile-only header失敗を含む凍結証拠を変更せず、別REVIEWとして保存する。

入口は `code/s61_nonsquare_P100_K7W_finite15_review_20261104_v63.py`。
新しい保存先を指定して再現する。native input/executableはignored work内に置き、
巨大係数や元polynomial証拠をREVIEWへ重複保存しない。

独立checkerの初回nativeは、最後の旧五cofactorのmodular和だけで停止した。
その前の全81 LCD、独立Laplace N202、N=HQ+R、新QQ content/Bézout/二除算、
旧五H除算のassertionは全て完了していた。原因はchecker固有のlow-level
modular preinverse callbackであり、実行source・入力hash・最終assertionを保存した。
完了済みの行列式を繰り返さず、元QQ係数をPython整数の通常の剰余と逆で再読込みし、
旧五の保存済みreduction全係数と既存Bézout和1を独立検査して閉じる。
新しい再現用sourceのcallbackも通常の整数積へ修正した。
数学的入力・制限・対象minorは変更していない。
