# actual honest T_nil：m/n/h の有限 source operator と完全128×972表示

2026-11-04。零次数 actual nilideal (N=\operatorname{Nil}(D_0)\) と reduced-normalizer
cokernel \(C'=\operatorname{coker}(D_0\to E_{4,\rm red})\) を用いる。
以下 \(S=R_0/(\ell)=k[[z_{A,0},\ldots,z_{A,6},r_A]]\)。r_AとDVR parameter ellを区別する。

本packetは actual coefficient functionsに対する **有限な形式函数operator** を構成する。
scalar係数や無限級数を展開せず、823680個の mX scalar積を列挙しない。
新しい source/gauge solve、jet、点、GB、nil-generator kernel solve は行わない。
operatorは derivativeと正整数Euler-degree divisionを含むので、polynomial DAGやpure
unit-Hensel systemだけと称しない。未指定の係数や非正則なdivisionを自由変数にしてもいない。

## 1. 全sourceと固定target

受理済み完全 radical comparisonは

\[
 R_0^{780}\xrightarrow{X}R_0^{44}\xrightarrow{V}R_0^3,
 \qquad H_1=N,\quad H_0=C'
\]

を与える。Xは全640 actual relationsと全140 SyzQradを含む。
全係数はsame physical R0の函数であり、ell0へのpullbackはsame actual F9 unit inverse
\((z7,r_A,s_A=r_A)\)で定義する。旧CB r/g、CA gAと全source函数を変更・再除算しない。

新targetは
\(\gamma=(\tau_A+11\tau_B/12,\tau_A\tau_B,\tau_B)\)。old targetからのconstant matrix Pは
既存 Koszul inputに保存され、\(V_\gamma=P^{-1}\bar V\)。
実 \(\bar C'\) の三sector annihilatorsは

\[
 (z7,r),\quad(z7,r),\quad(z7,r^2).
\]

完全minimal Koszul直和の最初の三写像を
\(d_1:S^{24}\to S^3,d_2:S^{84}\to S^{24},d_3:S^{168}\to S^{84}\) とする。

## 2. canonical m と全divisionのregularity

fが \((z7,r^n)\) に入るとき、zの最初のi個だけを残すprefixをf_iと書き、

\[
 f=\sum_{i=0}^6z_i\frac{f_{i+1}-f_i}{z_i}+r^n\frac{f_0}{r^n}
\]

を使う。prefix差のcoordinate divisionは恒等的にregularで、最後のr-power divisionは
complete cyclic cokernelにおけるmembershipからregularである。
旧B-sectorでは \(y_i=z_i+c_ir^2\) というtriangular automorphismを先に使い、最後に戻す。
原点で非単元のrやconductorをinverseとして使わない。

actual \(V_\gamma\) の全44列・三成分へこのoperatorを適用し、literal
\(m:S^{44}\to S^{24}\) を構成した。\(V_\gamma=d_1m\) は全source上の等式である。
そのregularityの根拠は accepted **全次数の** \(\bar C'\) のdiagonal presentationであって、
有限jetの零性の推測ではない。

## 3. source-bound n と自動的なunit24

旧 fullV56 packetはactual \(n_{\rm old}:S^{24}\to S^{56}\) を持ち、
\(\bar V_{\rm old}n_{\rm old}=d_{1,\rm old}\)。old-to-new source mapは

\[
 Z=\operatorname{diag}(\operatorname{Res}_4(A^t),I_8),\quad
 \bar VZ=\bar V_{\rm old}.
\]

eight actual normalized nil imagesのsource liftは、新radical source columns
20〜23および28〜31というliteral unit vectorsである。
actual old Cbar8の四つの独立像(n1,n1*tauA,n1*tauB,n2)を使い、24個の新canonical
relationのold Cbar residueを定数結合で表した。全8coordinatesの等号を固定QQ arithmeticで
確認する。差の実target函数は旧三sector idealに入るため、前節のcanonical divisionで
\(H:S^{24}\to S^{24}\) を作れる。

\[
 n_{\rm pre}=Zn_{\rm old}H+L_{44,8}B_{8,24},\qquad
 V_\gamma n_{\rm pre}=d_1.
\]

ここでLは上のeight literal source columns、Bは保存されたconstant lift matrixである。
旧sourceの全函数とnew radicalの実12→9 A-adapterを入力として保持する。

\(M=mn_{\rm pre}\) とすれば \(d_1(M-I)=0\)。complete minimal Koszul syzygiesの
全係数はmaximal idealに入るから \(M(0)=I_{24}\)。従って

\[
 n=n_{\rm pre}M^{-1},\quad mn=I_{24},\quad V_\gamma n=d_1.
\]

M inverseはこの証明でorigin-unitと分かったmatrix inverseであり、未知の非単元inverseではない。
scalar Mの数値やseriesを計算したという主張ではない。\(d_1M=d_1\) から
\(d_1M^{-1}=d_1\) となり、最後の等号も保持される。

## 4. actual780列ごとの完全Koszul lift

\(v=m\bar X_j\) は各sectorで \((z7,r^n)\)-cycleである。
**元の形式函数列そのもの**に次を適用し、有限normal formへ先に落とさない。

最後の成分v7をz7=0へ制限するとr^n*v7(0,r)=0。正則環ではr^nが非零因子なので
v7(0,r)=0となる。ordered Hadamardで \(v_7=\sum z_iH_i\) とし、
\(w_i=v_i+r^nH_i\) とすれば \(\sum z_iw_i=0\)。cross pair (i,7) のlift係数をH_iとする。
残る21pair(i,j)について

\[
 h_{ij}=(E_z+2)^{-1}(\partial_iw_j-\partial_jw_i),\qquad i<j<7
\]

を使う。係数のz-degree dに対しd+2で割るので、全て正整数divisionである。
Koszul differential \(\delta\) と形式的de Rham differentialについて
\(\delta\partial+\partial\delta=E_z+\text{exterior degree}\) が成り立ち、
\(d_2h(v)=v\) を全形式次数で与える。

このcontractionはS-linearではない。source coefficientを外へ出したり、Normal formの
representativeに作用させたりしない。780個の元列それぞれに適用した後、その780本を
S-linearに並べてmatrix \(h:S^{780}\to S^{84}\) とする。この順序なら

\[
 d_2h=m\bar X
\]

が正確に成り立つ。二つのn=1/2と全28pairを含む224個の固定polynomial fixturesで
actualtemplateの \(\delta h=v\) を照合し、S-linearityを仮定すると誤るfixtureも確認した。
一般の形式函数に対する根拠は上の恒等式であり、有限fixturesからの外挿ではない。

## 5. 完全128×972presentation

\[
 \Phi:S^{44}\oplus S^{84}\longrightarrow\ker V_\gamma,
 \qquad\Phi(s,b)=(I-nm)s+nd_2b.
\]

任意のcycleのm像はd1-cycleなので、この写像はsurjective。
\(\Phi(s,b)=0\) へmを掛ければd2b=0、従ってb=d3c、さらにs=ntとなる。
よってkernelは正確に(nS24,0)と(0,d3S168)で生成される。

actual X_jは \((X_j,h_j)\) へliftし、\(\Phi(X_j,h_j)=X_j\)。従って

\[
 T_{\rm nil}=\ker V_\gamma/\operatorname{im}\bar X
\]

の完全presentationは

\[
 S^{972}\xrightarrow{
 \begin{pmatrix}\bar X&n&0\\ h&0&d_3\end{pmatrix}}
 S^{128}\longrightarrow T_{\rm nil}\longrightarrow0.
\]

relation順序はactual780、section24、Koszul168。
\(\Phi\) を掛けるとfree44内では \((\bar X,0,0)\) であって、780列をfree44で
零としたのではない。actual source relationsをquotientしたときに零となる。
完全性は、任意の関係からactual780 liftを引き、残りを上のkernel \(\Phi\) で表すことで従う。

24unitの消去による104×948の存在は保証されるが、その20列complementを数値的に
選んだliteral matrixは本packetでは作っていない。保存された128×972のoperatorは
具体的で完全であり、未知のhや任意のsyzygy選択を入力にしていない。

## 6. 確定する長さの式と未決事項

比較free complexのH1はN、H0はC'。regular parameter ellによるbase changeでは
projective dimensionが1なので、H2が未知でも正確に

\[
 0\to N/\ell N\to T_{\rm nil}\to C'[\ell]\to0
\]

を得る。accepted actual NのDVR-finitenessと \(\mu_DN=L+3\)、
\(\dim C'[\ell]=2\) より

\[
 \boxed{\operatorname{length}T_{\rm nil}=L+5\in[16,22].}
\]

具体的なscalar係数、minimal生成数、exactL、nilidealの全generatorや三乗は計算していない。
NがIAIBにannihilateされるとも仮定しない。finite348 matrixへの非正当な置換もしていない。
次はこの完全operatorの正当な係数queryまたはminimalizationであり、typed表示の存在だけを
実値計算済みと扱わない。
