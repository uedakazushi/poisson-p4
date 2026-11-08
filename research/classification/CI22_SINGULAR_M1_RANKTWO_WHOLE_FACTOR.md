# singularM1・CI22全十源：元Jacobi ideal内のactual超平面因子

2026-11-08、PRIMARY E-FINITE-SURFACE。
**全Poisson点の実五次は超平面で割れる。実被約既約quinticのE枝は全空。**
全supportのparameter場合分けを解き直さず、既存の元350式／全126 actualPfaffianに
対する一回のQQ係数rowspace計算だけで証明する。

## 1. 全sourceとintrinsicな超平面

`NONNORMAL_CI22_RANK_TWO_SEPARATE_REVIEW.md` のmodel singularM1を対象にする。
canonical integral CI22は

\[
 q_1=x_0x_1+x_3^2,\qquad q_2=x_0x_2+x_4^2,\qquad S=V(q_1,q_2).
\]

全百五十tangency／intrinsic abelian-conormal条件をEuler商した全十sourceを
p0,…,p9とする。全kernel・free coordinate inverse・全PGL surface orbitの網羅性は
先行別作者receiptを入力とし、今回新しい小sourceへ制限していない。
π|Sの零divisor R、normal root、vertical条件は追加しない。

Sのsingular lineは{x0=x3=x4=0}。
その一般点でdq1=x1 dx0、dq2=x2 dx0なので、embedded Zariski tangent hyperplaneは
一定の H=(x0=0) である。
従って以下のcanonical H因子を全PGL標準形から戻しても、単なるchosen coordinateの条件ではない。

## 2. actualPfaffian modulo 元Jacobi係数の全恒等式

全sourceの元projective Jacobiの350係数をJ0,…,J349とし、実signed Euler PfaffianをFと書く。
次の空間次数四の多項式を置く。

\[
\begin{aligned}
G={}&3p_2p_6q_1^2\\
 &+(-p_0p_3+p_3p_6+6p_4p_8-6p_5p_6)q_1q_2\\
 &-3p_1p_5q_2^2.
\end{aligned}
\]

G自体は空間次数四、parameter次数二である。厳密に

\[
 F-x_0G=\sum_{j=0}^{349}H_j(x)J_j(p),\qquad H_j\in\mathbb Q[x]_5.    \tag{1}
\]

が成り立つ。右辺の各元row係数は同名reportの
`all126_Pf_original350_lifts` とmonomial順序から直接読める。
Fの各degree五係数についてQQ定数のoriginal-row functionalを保存したため、
全126の恒等式を独立に再構成出来る。

具体的には、保存された元14-dimensional quadratic rowspace basisを元350行へ
QQ線形結合で戻し、全126Pf係数を同じ保存basisで還元した。
全residualは旧producerの保存済みremainderと完全一致。
残る十空間monomialをまとめると正確にx0Gとなる。
新しいparameter Groebner basis、未知方程式、rank minor又はJacobi operatorではない。

式(1)はparameter環全体のpolynomial identityであり、p8,p5その他を反転しない。
正規形のgeneric係数又は一点でFをfactorして全源へ外挿したのでもない。
全Poisson閉点でJj=0なので

\[
                              F=x_0G.                  \tag{2}
\]

よってrank四実被約既約quinticは存在しない。
Gが零となる全補集合ではF=0、rank二以下であり、rank四sourceとして捨ててよい。
schemeの非被約厚み又は全Poisson idealのradical性は結論に不要である。

## 3. factor型への全点の戻し方とscope

Gは二quadric q1,q2のbinary quadraticである。
従って複素数体上の各rank四点ではGは二つのquadric pencil memberの積になる。
Sはintegral非退化なのでpencil memberのrankは三以上である。
actualDの全pointは **H+Q+Q'** 又は **H+2Q** という既存の最上位factor義務へ戻る。
double-Q型の全点既知被覆は既受理の全H+2Q theoremを使える。
distinct-Q型の全点が既知成分に入ると、本証明だけから追加主張しない。
その全source認識は既存の1+2+2 frontierに属し、ここで新しい曲面sourceの穴は作らない。

本稿が閉じるE-finite-surfaceの義務は、**actualDが被約既約quinticであるsingularM1 rank二枝の空性**。
normal CI22、regular22、nil32、cone又はprojected Veroneseをこの一modelと同一視しない。
新componentの増分は零、全E／全分類完了も主張しない。

## 4. 証明書

`jobs/ci22_singular_m1_ranktwo_whole_factor.json` のgate PASS後、
全350／126の保存配列に対して上の有限計算を行った。
**0.229637秒、RSS51,712KiB、全exact PASS。**
元14basisの全original-row lift、全126functional、全residual、全二次parameter monomial順序、
実行inline source、入力SHAは
`CLASSIFICATION_CI22_SINGULAR_M1_RANKTWO_WHOLE_FACTOR_20261108.json` に固定した。
個別caseの長いsupport cover、新Pfaffian profile、GB、minor又はjetを追加していない。
旧code・旧report・currentjob・主稿・global台帳・Gitは変更していない。
