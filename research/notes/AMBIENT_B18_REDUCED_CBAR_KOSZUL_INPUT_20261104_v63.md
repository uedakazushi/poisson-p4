# actual reduced Cbar の三 sector Koszul 入力と honest nil-comparison の source gate

2026-11-04。既存 cache の実値とその独立監査が与えた

\[
 \bar C'=k\gamma_1\oplus k\gamma_2\oplus k[[r_A]]/(r_A^2)\gamma_3,
\]

\[
 \gamma_1=\tau_A+\tfrac{11}{12}\tau_B,\quad
 \gamma_2=\tau_A\tau_B,\quad\gamma_3=\tau_B
\]

を、実 S=R0/(ell)=k[[zA0,...,zA6,rA]] 上の完全な minimal free resolution として保存する。
元の R0-coker C'、元 nilideal N、または ell-torsion の full resolution と同一視しない。

各 sector の regular sequence は

\[
 (z_0,\ldots,z_6,r),\qquad(z_0,\ldots,z_6,r),\qquad
 (z_0,\ldots,z_6,r^2).
\]

三つの Koszul complex の直和なので exact。全係数は maximal ideal に入り minimal である。
free ranks は

\[
 (3,24,84,168,210,168,84,24,3).
\]

同じ sector/wedge 順序で全3072個の signed polynomial entries を literal に保存し、
全 adjacent differential の積が零であることを固定 monomial arithmetic で確認する。
特に d1 は3×24、d2は24×84、d3は84×168。

## Target と比較行列の向き

元の target 順序は \((\tau_A,\tau_B,\tau_A\tau_B)\)。新 basis の列を並べた matrix は

\[
 P=\begin{pmatrix}1&0&0\\11/12&0&1\\0&1&0\end{pmatrix},\qquad\det P=-1.
\]

従って actual reduced-comparison 行列 Vrad の新 target 表示は \(P^{-1}\bar V_{\rm rad}\)。
次の factorization の向きは

\[
 P^{-1}\bar V_{\rm rad}=d_1m,\qquad
 \bar V_{\rm rad}n=P d_1.
\]

Pを source側に掛けたり、C/R0 restriction 内で転置したりしない。

## 次の source section の正確な入力

旧 actual Cbar resolution の typed m24×56、selected24 unit block M とそのinverseは
`FULLV56_COMPARISON_UNIT` に保存されている。そこから n_old56×24 と
\(\bar V_{\rm old}n_{\rm old}=d_{1,\rm old}\) が正確に得られる。
旧48normalizer columnsから新36radical columnsへの source mapは

\[
 Z=\operatorname{diag}(\operatorname{Res}_4(A^t),I_8),\qquad
 \bar V_{\rm rad}Z=\bar V_{\rm old}.
\]

八つの normalized nil image の lift は、新36列の indices20〜23および28〜31という
literal source unit vectorsである。actual Cbar8内で独立なのはそのうち
\(n_1,n_1\tau_A,n_1\tau_B,n_2\) の四本。各新 canonical relation の古い Cbar classを
この四本の定数結合で表し、その source liftを引く。残る exact target vectorは既知の
旧sector ideals \((z7,r^3),(z7+c r^2,r^3),(z7,r^2)\) に入るため、ordered coordinate
Hadamardと最後のregular r-power divisionで d1_old 係数を正確に作れる。
これは有限値から未知高次函数の零を仮定する操作ではなく、**全次数の既知cokernel**における
membershipによるregularityである。

新m24×44も同じ canonical divisionを新三sector idealsに適用して定まる。
section nを作れば \(d_1(mn-I)=0\)。minimal Koszul syzygyは全係数が maximal idealに入るので
\(mn\equiv I\pmod{\mathfrak m}\)。従ってunit inverseでnを正規化してmn=Iにできる。
この段階の具体的な m/n/quotient roots は本packetではまだ構成しない。

## 数値complementを選ばない完全な次のpresentation

mn=Iが source-boundに構成された場合、\(\ker\bar V\) を

\[
 \Phi:S^{44}\oplus S^{84}\to\ker\bar V,\qquad
 \Phi(s,b)=(I-nm)s+n d_2b
\]

で生成できる。その kernelは正確に \((nS^{24},0)+(0,d_3S^{168})\)。
各 actual780 relation X_jの mX_jはd1-cycleなので、完全Koszul contractionで
\(d_2h_j=mX_j\) と解けば \((X_j,h_j)\) が X_jを liftする。
従って m/n/hが構成された後は、honest T_nil は target128/source972で表示できる。
24unitを後で消せば target104/source948になる。complementのorigin数値を未確認のまま
選んだとはしない。

ellregularityとNilD0のDVR-finitenessは別の受理済みsource theoremにより閉じたが、
本packetはこの m/n/h780 constructionそのもの、元 nilidealの全生成元、または nilcubeを
完成したと主張しない。rA と DVR parameter ellを区別する。
