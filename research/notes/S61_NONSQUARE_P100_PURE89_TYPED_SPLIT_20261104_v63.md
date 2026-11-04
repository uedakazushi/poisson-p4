# P100 の固定九 target の全域分裂と pure89 同一 cokernel

2026-11-04。source-bound な有限 typed constructor。
新しい scalar 係数・normal form・rank を計算しない。下記の QQ Bézout、
CRT integer lift、norm15、adjugate、exact quotient の各 coefficient operator
は具体的な有限アルゴリズムとして保存するが、**係数値には展開していない**。
全曲線239や source kernel の自由基底は主張しない。

## 固定 projection と向き

既存 full \(\Gamma\) は319×98の**行** source 行列である。混乱を避け、
本記録では column map \(A=\Gamma^t:R^{319}\to R^{98}\) を使う。
\(R=T_{100}=\mathbb Q[z,1/(z(z+20)(z+60)(203z+2700)(z+12))]\)。

旧五個の9-minorの source binding を全て確認した。いずれも

\[
C=(0,1,2,3,4,5,6,9,25)
\]

という同一の九 target を使う。最後の列は original122 の30、parent239の
115であり、新 K7W 検査の追加列28とは異なる。\(E_C:R^9\to R^{98}\) を
この標準埋込み、\(\Pi=E_C^t\)、残る89列の埋込みを \(E_T\)、
\(\Pi_T=E_T^t\) とし、\(F=\Pi A\) と置く。

旧五個の unit-normalized minor の ideal は厳密に `(H15)`。
従って `D(H)` では F は全射である。全 H15 上では受理済み fullΓ98
証明から F の全射が従う。以下は、この二つを有限演算でつなぐ構成である。

## H 倍の右逆

既存15 sources の実319 sources への係数行列を \(J:R^{319}\to R^{15}\)
と書く。すなわち旧 \(M_{15}=J\Gamma\)。各 case の九行を、その**九列全体**
の保存済み unit LCD で清算する weighted selector を \(W_i:R^{15}\to R^9\)
とし、\(S_i=J^tW_i^t\), \(X_i=F S_i\) と置く。厳密な保存恒等式は

\[
\det X_i=U_i H Q_i,
\qquad U_i\in R^\times,\qquad Q_i\in\mathbb Q[z].
\]

U は保存された非零 QQ content と五つの既知 unit factor の冪だけである。
全五個の literal H-divisions と degree-preserved mod32003 cofactor Bézout
から \(\gcd_{\mathbb Q[z]}(Q_1,\ldots,Q_5)=1\) が厳密に受理されている。

\(\sum_i\beta_i Q_i=1\) の polynomial coefficients は未計算である。
constructor は固定した i の順で monic QQ[z] extended Euclid を行い、
以前の重みを更新するという canonical finite operator として \(\beta_i\)
を指定する。この有限演算を実行したと装わない。

\[
R_H=\sum_i\frac{\beta_i}{U_i}\,S_i\operatorname{Adj}(X_i)
\quad\Longrightarrow\quad F R_H=H I_9.
\]

H や Q_i を regular unit として逆にしていない。

## 有限候補上の adj98 と norm による接合

七因子で受理された selected98 source IDs は、全域の **実 A の source
columns** を固定する。それぞれの埋込みを \(S_j:R^{98}\to R^{319}\)、
\(B_j=A S_j\)、\(D_j=\det B_j\) とする。B_j は実 global Γ から得る
lazy matrix であり、有限体の係数を適当に QQ へ lift した行列ではない。

\[
R_j=S_j\operatorname{Adj}(B_j)E_C,
\qquad A R_j=D_jE_C,
\qquad F R_j=D_jI_9.
\]

Hbar の七因子に対する CRT idempotents は旧 inventory には係数として
保存されていない。ここでは \(h_j=\bar H/f_j\) とその canonical Fp inverse
を使い、`reduce(h_j*inverse(h_j mod f_j),Hbar)` を求め、各係数を
0..32002に lift する有限 operator として \(e_j\in\mathbb Z[z]\) を指定する。
これも未実行であり、保存済み idempotent との偽った記述はしない。

\[
G=\sum_j e_jD_j,\qquad R_{\mathrm{fin}}=\sum_j e_j R_j.
\]

\(\bar G\bmod f_j=\bar D_j\bmod f_j\ne0\) は、独立に確認された全七個の
非零98-minorから従う。従って \(\bar G\) は全 \(\bar S\) の単元である。
\(S=\mathbb Q[z]/H\) の基底 \(1,z,\ldots,z^{14}\) で G の乗算を表す
15×15行列を \(M_G\) とする。係数の curve denominators は H と互いに素
なので、この行列も canonical quotient-ring 演算で構成できる。

\[
n=\det M_G\in\mathbb Q^\times,\qquad
\mathbf b=n^{-1}\operatorname{Adj}(M_G)\mathbf e_0.
\]

norm の reduction は CRT で全成分が可逆なため非零であり、n は非零 QQ
constant（実際には32003-unit）。これは z に依存する未知 factor の逆では
ない。\(\mathbf b\) から次数15未満の polynomial b を得て、\(bG=1\bmod H\)。

`1-bG=N/D` と既知の曲線 unit のみからなる共通分母で書けば、H と D は
互いに素であり H は N を文字どおり割る。従って

\[
a=\operatorname{ExactQuotient}_H(1-bG)\in R,
\qquad aH+bG=1.
\]

この exact quotient は numerator clearing と QQ[z] の exact polynomial
division で指定する。H を逆にして0/0から regularity を推測する操作ではない。

\[
\boxed{R_9=aR_H+bR_{\mathrm{fin}},\qquad F R_9=I_9.}
\]

各 adj9/adj98/adj15 は universal polynomial operator のまま保持する。
個々の未知 D_j や H を全域で invert する必要はない。

## 全319 sources を残す pure89 presentation

\(V=\Pi_T A R_9\) とする。C,T の順に target を並べれば、

\[
[A R_9,E_T]=\begin{pmatrix}I_9&0\\ V&I_{89}\end{pmatrix},
\quad \det=1.
\]

元の target 順では既知の並べ替え符号が加わるだけである。
source を310個の自由基底に変換したとは主張しない。全319 generatorsを
そのまま残し、

\[
P=I_{98}-A R_9\Pi,\qquad
T_{89}=\Pi_T P A:R^{319}\longrightarrow R^{89}
\]

を使う。P は \(\ker\Pi\) への projection、\(T_{89}R_9=0\)。
target map \(\Pi_TP\) と逆向きの \(E_T\) が cokernel 上で互いに逆なので

\[
\operatorname{coker}(A)\simeq\operatorname{coker}(T_{89}).
\]

行 source 表示では319×89の同一 cokernel presentation となる。
元の χ0、全 cubic pure tails、全33 quartic source は残る。98 scalar
targets の展開、89行列の評価、source kernel310の自由性、残る89列の全射は
一切仮定しない。従って全曲線239は依然として未解決である。

## 保存形式

`code/s61_nonsquare_P100_pure89_typed_split_20261104_v63.py` は、凍結 gate の
hash/source/target/indexと全 operator shape を照合し、有限 operator graphを
`reports/S61_NONSQUARE_P100_PURE89_TYPED_SPLIT_20261104_v63/operator.json`
に保存する。新しい parameter polynomial/RF coefficient terms は0。
canonical coefficient operators の未実行状態を report/proof に明記する。
独立レビューは別 gate とする。
