# nil22 image double_B全六磁気源のCPLUS_EXT帰属：独立検証

2026年11月7日、主義務E。**判定：全源についてPASS。**
`nil22_image_double_B` の全六磁気係数を持つ源は、既知24次元成分 `R4_CPLUS_EXT` に含まれる。
新成分は加えない。この一源の帰属を閉じ、他の五つの未認識源は残す。

独立実装は `code/classification_nil22_double_b_cplus_limit_review.py`、
証明書は `reports/CLASSIFICATION_NIL22_DOUBLE_B_CPLUS_LIMIT_REVIEW_20261107.json`、
有限契約は `jobs/nil22_double_b_cplus_limit_review.json`。
producer・CAS・repository arithmeticをimportしないstandard-library Fraction計算である。
契約gate PASS後、全恒等式は約0.0074秒、最大RSS12,032 KiBで確認した。

## 1. 全源を保つ実曲線とendpoint

座標を(h,q₁,q₂,x,y)、εを曲線parameterとし、
R=m₀h²+m₁hx+m₂hy+m₃x²+m₄xy+m₅y²とする。
t=ε−2、δ=4tと置く。h-Casimir二ベクトルの全非零括弧の共通分子N=δPは

\[
\begin{aligned}
 N^{12}&=4tR-(ty+4(t+4)h)q_1,\\
 N^{13}&=-tx^2+8h^2,&N^{14}&=-2txy-8(t+4)hx,\\
 N^{23}&=-txy-4(t+4)hx,&N^{24}&=32th^2-2ty^2+32(t+2)x^2.
\end{aligned}
\]

N³⁴とhとの括弧は零である。
全m₀,…,m₅を独立変数のまま保持し、producerの全斉次・affine表示と直接一致することを確認した。
凍結した元の完全源 `CLASSIFICATION_NIL22_IMAGE_REMAINING_SOURCES_20261107.json` の
case1、`nil22_image_double_B` について、全十括弧で

\[
 N|_{\varepsilon=0}=-8P_{\mathrm{original}}
\]

となる。δ(0)=−8なのでP=N/δのendpointは元の全源そのものである。
全四strict affine Jacobiと、全五四添字・全七十単項式の350projective Jacobiを
別に構成し、全て恒等零と確認した。
六磁気係数の一つも反転又は特殊化していない。

## 2. 正しい実Pfaffianとgeneric配置

\[
 Q_A=hy-4h^2+\varepsilon x^2,\qquad Q_B=hy+4h^2+2x^2.
\]

全五つの符号付きPfaffianからEuler収縮を直接計算し、

\[
 D_N=16(\varepsilon-2)hQ_AQ_B,\qquad
 D_P=\frac{hQ_AQ_B}{\varepsilon-2}
\]

と確認した。三変数quadratic matrixの行列式はそれぞれ−ε/4、−1/2である。
従ってε≠0,2で二conicはrank三・既約・相異なり、H=(h=0)とも異なる。
Hへの制限はεx²、2x²なので、Hは共通点[0:0:1]での共通接線となる。

実可逆線形変更(h,x,y)=(X,Y,−2Z−4X)によって

\[
 Q_B=-2(XZ-Y^2),\quad
 Q_A=-\varepsilon\bigl((2/\varepsilon)XZ-Y^2+(8/\varepsilon)X^2\bigr).
\]

これはshared_tangent_genericで、u=2/ε≠0,1、d=8/ε≠0である。
ε=4のu=1/2も、既受理の別kernel chartに含まれ、除外しない。

## 3. 全六磁気定理への適用は因子形だけの認識ではない

Nの二つのbase anchorをN_U=(N¹³,N¹⁴)、N_V=(N²³,N²⁴)とする。
独立計算では、単なるfactor labelに加えて全log-tangencyを検査した。

\[
\begin{array}{c|cc}
 &Q_A&Q_B\\\hline
 N_U(Q)&-2txQ_A&-2txQ_B\\
 N_V(Q)&-2t(y+4h)Q_A&-2t(y-4h)Q_B
\end{array}
\]

またdet(N_U,N_V)=−16tQ_AQ_Bなので、Eulerと二anchorのSaito determinantは
非零scalar倍のhQ_AQ_Bである。
したがってgeneric配置の完全free frameに属し、chosen pairだけを全支持と仮定していない。

上のbase変更後の磁気二次式の六係数は

\[
 (m_0-4m_2+16m_5,\ m_1-4m_4,\ -2m_2+16m_5,\ m_3,\ -2m_4,\ 4m_5)
\]

となる。元の六係数からの変換行列式は16で、全磁気A⁶を保つ。
残る一定fiber frame変更もfree-frame bundleの全frameとして、既受理の全源定理に含まれる。

`TWO_QUADRICS_PLANE_CONIC_CONFIGURATIONS.md` §§1,4–5の
shared_tangent_generic **全六磁気・全free-frame源**がCPLUS_EXTへ入るという幾何定理を適用する。
そこでは既知唯一成分点を含む既約frame bundleの像閉包を用いており、
因子配置一例又はgeneric次元だけから帰属させたものではない。
本reviewはその既受理定理を入力とし、今回の全curveが全仮定を満たすことを確認した。

## 4. 閉性と補集合

任意に固定した六磁気係数について、P_εはε≠2の実正則Poisson族である。
ε≠0,2の全fiberは上の全源定理からCPLUS_EXTに入り、
CPLUS_EXTの閉性によりε=0の元の源も入る。
磁気係数に条件が無いため、その全支持とPGL飽和閉包が含まれる。

ε=2はこの曲線の極であり除外するが、endpoint閉性には不要である。
εの反転はgeneric conicの識別だけに使い、ε=0を消していない。
他のnil22 image源、主分類又は全Eの完了は主張しない。
