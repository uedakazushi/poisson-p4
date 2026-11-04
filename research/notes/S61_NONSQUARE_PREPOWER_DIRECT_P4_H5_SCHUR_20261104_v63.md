# h22 を反転しない実 direct P4＋H5 の6×13 Schur表示

2026-11-04。新 scalar polynomial/RF coefficient terms は0。
受理済み pre-POWER の実15×22入力に対する source-bound typed operator。
残る6×13の値や階数は計算しない。

## 一貫した行・列と単位領域

実 P4 は current15 source rows `[0,1,3,6]`、すなわち
`[Q9,G18,G39,G68]`、original32 targets `[22,11,13,25]` である。
独立監査済みの polynomial identity から

\[
D_4=\det P_4=a^4 C^3 R_{434}.
\]

これは h22 消去後の3×3 determinant ではなく、その**前**の実4×4 block。
逆は `Adj4(P4)/D4` を直接使い、h22 を一度も反転しない。

H5 の compact14 row indices `[1,3,6,12,9]` は current15 では
`[2,4,7,13,10]` であり、P4の四行とは重複しない。
heavy columns は compact5順 `[0,4,1,3,2]`、original32順では
`[27,31,28,30,29]` である。

既存70恒等式は `Hactual=diag(rho) N diag(1,1,1,1,gamma)`。
P4 determinant を文字どおり上の D4 に保つため、**H5に選んだ五source行だけ**
rhoで割り、最後のheavy target31をgammaで割る。他の十source行は未変更。
rho,gamma は既知FULLU unitsである。これは全14行を再正規化する操作ではない。

以後の領域は正確に

\[
\mathrm{FULLU}\cap D(R_{434}\Delta_5).
\]

R434とDelta5の零点側は残す。h22=0を追加で除外しない。
R434は未知因子であり、全FULLU単位や欠損 locus とはしない。

## H5 は既存の2＋3 cofactor blockで逆にする

選択された normalized H5 は、実 source index/coefficients から

\[
H_5=\begin{pmatrix}P_2&0\\ C_{32}&B_3\end{pmatrix},\quad
\delta_2=\det P_2=-\tfrac65 Bp^2,\quad
\delta_3=\det B_3=\Delta_{407},\quad
\Delta_5=\delta_2\delta_3
\]

である。上右六成分は実 literal zero。B3 の row/col順は既存rank証明書の
selected rows7 `[0,5,2]` / cols3 `[0,2,1]` に一致する。
1752項の既存Delta5を再展開しない。

\(A_2=\operatorname{Adj}_2(P_2)\), \(A_3=\operatorname{Adj}_3(B_3)\) とすれば

\[
H_5^{-1}=\frac1{\Delta_5}
\begin{pmatrix}\delta_3 A_2&0\\-A_3C_{32}A_2&\delta_2 A_3\end{pmatrix}.
\]

両側積は2×2、3×3の既存 cofactor identity だけで Delta5 I5 となる。
新しい5×5随伴計算をしない。2/3のuniversal cofactor operatorsは未展開であり、
保存されていないnumeric adj3を計算済みと呼ばない。明示的な1/pや1/Delta407も
追加せず、inverse nodeはDelta5だけである。

## 九 pivot と実6×13

並べ替えた pivot source順は
`[0,1,3,6,2,4,7,13,10]`。
target順は original32 の
`[22,11,13,25,27,31,28,30,29]`。
対応する9×9 blockは

\[
K=\begin{pmatrix}P_4&0\\ C_{54}&H_5\end{pmatrix},\qquad
\det K=D_4\Delta_5,
\]

\[
K^{-1}=\begin{pmatrix}
P_4^{-1}&0\\-H_5^{-1}C_{54}P_4^{-1}&H_5^{-1}
\end{pmatrix}.
\]

残る source rows は current15 の `[5,8,9,11,12,14]`。
残る13 seed targets は original32 の
`[1,2,3,4,6,7,8,9,14,15,17,18,19]`。
Q9の、選択h9九係数以外の**八つ全てのpure tail**がこの13列に含まれ、
B block に実 root のまま保持される。Schur後の列はseed列の単なる部分行列
ではなく、pivot列を加えた実 target transform の結果である。

full matrix を `[[K,B],[C,D]]` と書き、

\[
L=\begin{pmatrix}K^{-1}&0\\-CK^{-1}&I_6\end{pmatrix},\qquad
U=\begin{pmatrix}I_9&-K^{-1}B\\0&I_{13}\end{pmatrix}
\]

とすれば、

\[
L\begin{pmatrix}K&B\\C&D\end{pmatrix}U
  =\operatorname{diag}(I_9,\ D-CK^{-1}B).
\]

L/Uの逆、元のhybrid行単位/列単位と並べ替えを含む full source/target inverse
も保存する。従って元のactual15×22と新actual6×13は、指定領域で同一cokernel。
presentationの単位九列を消しており、未知projective kernelの自由基底を
仮定していない。rank は `9+rank(residual)` だが、residual rank6は未証明。
full cubic57やquartic239もここから結論しない。

## 元25 sourceへの係数輸送

六行の残存source covectorを current15 へ明示し、既存 source16/source17
unit matrices、各削除行への零の挿入、八Gauss source stepの順序付き積を通して
元25 sourceへ戻す。prefixは **forward** row operatorの合成である。
`rowr:=rowr/p; rowt:=rowt-c_t rowr` のstepを逆順にcoefficientへ作用させる
ことと同じであり、sourceをtransposeしたり無関係な逆行列にしたりしない。
これらpは全て既知FULLU units。h22/Delta/Hadj/POWER207の係数は輸送に使わない。

保存コード:
`code/s61_nonsquare_prePOWER_direct_P4_H5_schur_20261104_v63.py`。
同名 report の operator.json は全 matrix shapes、source/target permutations、
unit inverses、残存6sourceの元25への輸送、全tail参照を固定する。
実 parameter products、scalar minor、NF、rank、source点、prime、GBを新たに
計算しない。独立レビューは別gateとする。
