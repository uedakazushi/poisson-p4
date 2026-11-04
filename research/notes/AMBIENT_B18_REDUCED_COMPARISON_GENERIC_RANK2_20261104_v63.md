# actual reduced-normalizer cokernel の generic DVR rank2

2026-11-04。\(\psi:D_0\to E_4\)、\(C=\operatorname{coker}\psi\) とし、

\[
 C'=\operatorname{coker}(D_0\longrightarrow E_{4,\rm red})
\]

を同じ physical R0 上で定義する。accepted comparison の DVR は
\(D=k[[\ell]],\ell=5(s_A-r_A)\) である。本稿の結論は

\[
 \boxed{\operatorname{rank}_D C'=2},\qquad
 C'\simeq D^2\oplus D/(\ell^{m_1})\oplus D/(\ell^{m_2}),\quad m_i\ge1
\]

という非標準的な D-module 分解である。二つの指数や torsion の長さは未計算。
元 nilideal が DVR-finite であることや、その生成数をこの分解から推論しない。

## 1. 既存 C4 の same-source inventory を quotient に移す

\(E_4\to E_{4,\rm red}\) は全 nilideal による quotient なので、C' は C の quotient。
従って C を annihilate する \(I_AI_B\) は C' も annihilate する。
既受理の intrinsic conductor と exact same-base joint ideal theorem は、この零次数底に
CA-critical、CA-self Pa、CB-critical、CB-self Qb の四つの generic support だけを与える。
C' のために新しい支持枝の尽くしを仮定しない。C' が C の quotient であることによる
支持の包含を使う。

各 support curve は smooth で、同じ source の \(\ell\) は order1、generic reduced field
の D-degree は1である。これは sourceC4 theorem とその独立監査が、actual CA/CB
conductor embedding と proper-fiber maps を束縛して証明したものである。

## 2. full source から zero comparison への輸送

共通 actual G19 action-product は、両 normalization inclusion とその tensor map、cokernel
に同時に作用する。regular formal group factor は flat であり、剰余次元20の source 模型の
generic 長さを横断底へ移しても変えない。別の abstract normal form の座標を physical R0 と
同一視したり、非 flat な断面を勝手に tensor したりしていない。

generic 零次数 conductor prime では、相手の conductor が単元で、当該側の positive J_i0
も単元となる。実 \(I_i=I_{i0}+wJ_{i0}\) から w は当該 full conductor に入り、full C は w に
殺される。源の coker 生成元 \(\tau_i,\tau_A\tau_B\) は character0 なので、この full coker は
degree0 に集中する。これが sourceC4 packet で既に受理された厳密輸送である。

full reduced-target coker はその graded quotient だから、同じ集中性を保つ。
非負gradingでは nilradical の degree0部分が degree0 ring の nilradical と正確に一致するため、
\((E_{\rm full,red})_0=(E_4)_{\rm red}\)。degree0を取る操作は graded module に対して exact。
従って下の full source quotient-map 計算は、この同じ C' の generic localization を計算する。
単に w に殺される module を degree0と推測した議論ではない。

## 3. actual四写像の reduced-target cokernel

* **CA-critical。** accepted actual contact は
  \(D=G[[l_7,v,x]]/x(l_7,x,v^2)\)、
  \(E=G[[l_7,v,t]]/(l_7t,t^2-vt,t^3)\)、\(x\mapsto t^2\) を与える。
  両 reduced quotient は同じ regular \(G[[l_7,v]]\) で、写像は恒等。
  従って C' の寄与は0。
* **Pa。** accepted whole proper normalization fiber と全 actual scheme intersection は、
  reduced normalized tensor を \(\mathcal O_N\oplus\mathcal O_{S_P}\) と同定する。
  元 D の reduced quotient は \(\mathcal O_N\)、第二成分はその closed reduced S_P への制限。
  よって diagonal graph の cokernel は \(\mathcal O_{S_P}\) で、generic長さ1。
* **CB-critical。** accepted whole proper-fiber/global boundary theorem により actual CA の外。
  元比較も quotient 比較も零で、寄与0。
* **Qb。** accepted asymmetric whole source map は、CAのregular環 B と
  \(J_Z=(z_7,f^2,\rho^2)\) を使い
  \(D=B/(fJ_Z)\to E=B/(f)\oplus B/J_Z\) という二つの quotient map である。
  reduced target は \(B/(f)\oplus B/(z_7,f,\rho)\)。後者への写像は前者への quotient に因子化する。
  従って reduced diagonal graph の cokernel は \(B/(z_7,f,\rho)\)、generic長さ1。

ここでは Qb の whole nil module 長さ4や元 coker 長さ2を、reduced coker の長さ1と混同しない。
Pa/Qbはラベル交換による推測ではなく、各々独立監査された実 source map を使う。
\(A\to A\oplus A/I\) の graph の cokernel が A/I であることは、
\((a,b)\mapsto b-a\bmod I\) で直接分かる。

四つの寄与と各 curve の D-degree1より

\[
 \operatorname{rank}_D C'=0+1+0+1=2.
\]

## 4. 原点 fiber と torsion の生成数

新しい actual B11 値と literal Cbar sectors から、別 packet は

\[
 C'/\ell C'\simeq k\oplus k\oplus k[[r_A]]/(r_A^2),
 \qquad \operatorname{length}(C'/\ell C')=4
\]

を正確に与える。この (r_A\) は \(\ell\) を切った後の別の係数座標であり、DVR parameter
ではない。C' は C の quotient なので D 上有限である。DVR 構造定理から
\(\mu_D C'=4\)、free rank2を引いて torsion 最小生成数は2となる。

\(m_1,m_2\) は未知であり、有限 \(\ell=0\) 値から1などと決めない。
元 D0 の nilideal、K/JK の honest presentation、nil三乗・正確な nil 指数、全域の支持分類は
本稿の結論ではない。
