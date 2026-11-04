# 特殊比較 module T：正確な完全列と小さい source 比較手順

2026-11-04、v63。方法・source preflight であり、実 (T\) の大行列や長さを計算しない。

\(D=\mathbb Q[[\ell]]\)、実正則底 (R=R_0\)、\(\bar R=R/(\ell)\) とする。
既受理の actual (UP=QB\) と complete normalizer syzygy source を用い、
normalizer の12生成元を regular8＋零4へ unit 変換し、四成分 restriction で16零列を除く。
その後 actual unit1 の source/target を消すと、比較の低次数 complex は

\[
 R^{752}\xrightarrow{X}R^{40}\xrightarrow{V}R^3.
\]

元の符号は (V=(Q_{\rm red},-U)\) とし、(X\) に
\((B_{\rm red},P)\) の640列と \((\operatorname{Syz}Q_{\rm red},0)\) の112列を入れる。
この source 変換は既存の unit operator から定義できるが、本 packet は実 (X,V\) の scalar 行列を新たに保存したものではない。

\[
 J=\operatorname{im}V,\quad L=\ker V,\quad K=L/\operatorname{im}X,\quad C=\operatorname{coker}V
\]

と置く。ここで (K\) は元の actual 比較 kernel である。

## 正確な特殊化

(J\subset R^3\) なので (J\) 上の \(\ell\) は非零因子である。
二つの短完全列を tensor し、

\[
 0\longrightarrow L/\ell L\longrightarrow\ker\bar V
 \longrightarrow C[\ell]\longrightarrow0
\]

を得る。actual specialized (\bar X\) の像は左端の中にあり、その像で割れば

\[
 \boxed{0\longrightarrow K/\ell K\longrightarrow
 T\longrightarrow C[\ell]\longrightarrow0},\qquad
 T=\ker\bar V/\operatorname{im}\bar X.
\]

左端の同一視は cokernel の右完全性
\((L/\ell L)/\operatorname{im}\bar X=K/\ell K\) による。
特殊化 kernel を (K/\ell K\) と同一視しない。完全 Syz(P\) はこの完全列に不要である。

実 (C\) は有限 (D\)-module であり、既受理の
\(\mu_D C=8\)、\(\operatorname{rank}_D C=4\) から
\(\dim_{\mathbb Q}C[\ell]=4\) が従う。
従って、将来実 (T\) の長さが求まれば

\[
 \boxed{\mu_D K=\operatorname{length}_{\mathbb Q}T-4}.
\]

これは (C\) の torsion の全長が4という主張ではない。
また (K\) の一般 rank を決めない。
有限 (D\)-module (K\) の torsion が非零であることを判定するには、
この最小生成数を別途確定した \(\operatorname{rank}_D K\) と比較する必要がある。

## 実 Cbar 分解を使う完全な小比較

別 packet `CBAR_LITERAL_RESOLUTION` は、actual target
\((\tau_A,\tau_B,\tau_A\tau_B)\) を持つ実 \(\bar C\) の最小自由分解

\[
 \cdots\to\bar R^{168}\xrightarrow{d_3}\bar R^{84}
 \xrightarrow{d_2}\bar R^{24}\xrightarrow{d_1}\bar R^3\to\bar C\to0
\]

を全3072多項式成分で構成した。
八 basis の三 target への持上げは
\(e_A,re_A,r^2e_A,e_B,re_B,r^2e_B,e_{AB},re_{AB}\) と実 source で一致する。
一般 DVR 上の (C\) の分解を先に作る必要はない。

次に actual 関数係数を持つ比較

\[
 m:\bar R^{40}\to\bar R^{24},\quad
 n:\bar R^{24}\to\bar R^{40},\quad
 h:\bar R^{24}\to\bar R^{84}
\]

を構成し、

\[
 d_1m=\bar V,\quad \bar Vn=d_1,\quad d_2h=I-mn
\]

を確認する必要がある。そのとき

\[
 S=[nd_2,I-nm]:\bar R^{124}\twoheadrightarrow\ker\bar V
\]

は完全な生成写像であり、その完全関係は

\[
 \begin{pmatrix}d_3&-h\\0&n\end{pmatrix}:\bar R^{192}\to\bar R^{124}
\]

である。実際、(S(w,x)=0\) に対し (y=mx-d_2w\) とすると
\(x=ny\)、\(d_2(w+hy)=0\) となり、\(\ker d_2=\operatorname{im}d_3\) を使えばよい。

actual 各 (v=\bar X_j\) には \(d_2w=mv\) を満たす実持上げを選び、
\(S(w,v)=v\) とする。従って (T\) は124生成・192＋752＝944関係の
**完全な cokernel**として構成できる。単なる同値な kernel criterion ではなく、
必要な比較行列と cycle lift の形が固定されている。
ただし本 packet はその (m,n,h,w\) の係数を計算したものではない。

## 正しい有限底と残る gate

\(I\subset\bar R\) を actual (I_AI_B\) の像とする。
既受理の normalizer conductor 論から (IK=IC=0\) であり、上の完全列により

\[
 I^2T=0.
\]

(IT=0\) は従わない。実 (H_{20}=\bar R/I\) は \(\mathfrak m^5=0\) なので
\(\mathfrak m^{10}T=0\) も従うが、(H_{20}\) 上の行列の kernel を (T\) としてはいけない。
正しい有限底は \(\bar R/I^2\) である。
別担当の source-bound348環 packet を用いる場合、その source pair 正規化と独立 gate を保持する。

完全な (T\) cokernel ができてからこの有限底へ tensor する操作は右完全である。
保守的な未削減サイズは348×124＝43,152個の有理生成元、348×944＝328,512列なので、
このままの密行列は承認していない。適合する basis・constant cancellation・factored sparse operator を先に使う。

残る実 gate は、(X,V\) の source unit adapter、実 (m,n,h\) と752 cycle lift、
そして348環への全係数関数の正確な制限である。
元 (r/g\) が特殊化後に0/0となっても、新しく割り直したり零に置いたりしない。
各係数を実正則関数の restriction operator として保持する。

本 packet は (T\) の長さ、\(\mu_DK\)、\(\operatorname{rank}_DK\)、socle、元の深さを計算していない。
