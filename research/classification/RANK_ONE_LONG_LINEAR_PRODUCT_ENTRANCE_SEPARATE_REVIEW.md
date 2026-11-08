# rank一 long M1／M11：全25次元像と marked modular-rank 障害の独立純監査

2026-11-08、別作者 quartic_sources、PRIMARY E_multiple_hyperplane。**純全源入口 PASS**。
対象は `RANK_ONE_LONG_LINEAR_AND_PRODUCT_ENTRANCE.md`
（SHA256 `fcb2f23dde160b8d7a346acc42b6881f8c682230ba1c21f82fba963a040109ab`）。
新しい点、operator、Jacobi 行列、rank／minor、kernel、profile、jet 又は CAS は実行していない。
受理するのは二つの相異なる既約25次元 Poisson 像閉包、全既知からの分離まで。
両像の全空間成分性と未知の高次元親は、この入口だけでは未解決である。

## 1. 全 source の正確な次元と actual 群

旧全長鎖式(19)、\(X=\partial_x+x\partial_y\)、\(g=0\)、\((U+1)H=0\)、全六係数の \(F\) を使う。
\(H,U,F\) は tensor から一意に復元される。H の定数差も同じ固有値式で零となる。

M1 の H は任意の非定数 affine linear 式なので基底三次元。
H=z の actual affine frame で全解は
\(U=-z\partial_z+(v_0+v_1z+v_2w)\partial_w\)、三自由係数。
従って全基底上で非空 affine fiber の rank は三、source は既約 \(3+3+6=12\) 次元である。
Jordan、零固有方向の translation、全 F を捨てていない。

M11 の H は非平行二 affine line の積。
二本の affine line 各二次元と非零全体係数一で、H の基底は既約五次元。
線形部分が独立なので交点は一意。中心／frame で \(H=tzw\) にすると、
factor 別 semi-invariance から全 U は
\(\lambda z\partial_z+(-1-\lambda)w\partial_w\)、\(\lambda\) 全値。
よって source は既約 \(5+1+6=12\) 次元。
scalar／zero-eigenvalue／二 factor 交換の有限重複も同じ source に含む。
平行な相異二 affine factor は、各 factor を保つ affine U がその両方へ同じ非零 weight を与えられず、
\(U(H)=-H\) と両立しない。全有限 routing の M11 定義と一致する。

全 StabX 十次元は \(r'=r+\Delta(q)\)、\(q'=Bq+t\) として
\(U'=B_*U,H'=\det(B)H,F'=F-U\Delta-(1+d)\Delta\) を与える。
両 source は全域でこの actual 変換と逆により保たれる。
\((x,y,q)\mapsto(\lambda x,\lambda^2y,q)\) と tensor scalar \(1/\lambda\) による
\(H'=\lambda H,F'=\lambda^2F,U'=U\) も全域である。
作用写像の fiber は係数逆と divergence により正確に StabX の一 coset。
従って affine 飽和は22次元、scalar を除き H を動かす四次元を戻して25次元。

M1 の actual D は \(h^4L\)、M11 は \(h^3L_1L_2\)。
一般点で h は唯一の重複 hyperplane なので marked H を忘れる fiber は一点である。
両像閉包 \(Z_1,Z_{11}\) は正確に既約25次元となる。
標本の orbit 又は tangent rank は、この下界に使用していない。

## 2. Canonical marked modular matrix の全正則性

五次元線形空間の Euler 場を E、二次 skew lift の column divergence を
\(\delta(P)^j=\sum_i\partial_iP^{ij}\) とする。
\(\operatorname{tr}\delta(P)=0\)、\(\delta(E\wedge A)=5A-\operatorname{tr}(A)E\) なので
\[
 \overline\Pi=P-E\wedge\delta(P)/5
\]
は divergence 零であり、Euler gauge の選択に依存しない。
GL5 の constant-Jacobian 作用に equivariant である。

Poisson hyperplane \(\ell=0\) に対し **column-sharp**
\(Y=\overline\Pi\,\nabla\ell\) を取る。
\(Y(\ell)=0\)、projective tangency より \(Y\bmod\ell\) は Euler に比例する。
\(\ell=x_0\) に線形変更すると、次数二の全成分と \(x_1,\ldots,x_4\) の gcd 一から、
その比例係数は H 上の大域一次 section であり、一次式へ lift 出来る。
従って全域で
\[
 Y=\ell A-A(\ell)E
\]
と書ける。\(A\) は一次場。A の scalar-Euler の曖昧さは下式で消える。
divergence 零より
\[
 0=\operatorname{div}Y=\ell\operatorname{tr}A-5A(\ell).
\]
ゆえに \(B=A-(\operatorname{tr}A)E/5\) は trace零、\(B(\ell)=0\)、\(Y=\ell B\)。
\(\overline\Pi+E\wedge B\) が \(\ell\)-Casimir lift で、その divergence は
\[
                   N_\ell=5B=5\overline\Pi\nabla\ell/\ell.
\]
\(\ell=h\) の元 h-Casimir lift に戻せば、その元 divergence と正確に一致する。
row-sharp では符号が逆なので、以下の raw row cofactor 計算と混ぜない。

dual P4 の \(\ell_j\ne0\) chart で \(\ell_j=1\) とし、\(\ell\) を一つの座標に取る。
これは chart 上で可逆・regular な線形変更であり、\(Y\) の \(\ell\) による monic exact division は
係数に regular に依存する。従って \(N_\ell\) は全 marked Poisson hyperplane incidence 上で regular。
\(\ell\) の尺度に依存せず、tensor の尺度では同じ尺度が掛かるだけなので、
その matrix rank は projective tensor に対しても well-defined である。

actual D 非零の rank四開上では、\(\ell^m\mid D\)、\(\ell\) Poisson、
\(\operatorname{rank}N_\ell\le1\) は marked incidence の閉条件。
\(\ell\) が proper な dual projective space を動くため、存在条件も rank四開内で閉である。
Pf零への落下や発散する coordinate normalization から label を選んでいない。

## 3. CA の全 magnetic source と triple label

rawCA 実装を import／実行せず、保存された定義を読解した。
元全 source は
\[
 P_{xy}=f(x,h),\quad P_{zw}=g_0h^2,\quad
 P_{yz}=-H_{3,w},\quad P_{yw}=H_{3,z},\quad P_{xz}=P_{xw}=0.
\]
全十六 cubic monomial の任意係数を保つ。
mixed partial が相殺して \(\delta P=f_x\partial_y\)。
homogeneous matrix は出力が y の一行だけなので rank 一以下である。
generic actual D は \(3H+H_1+H_2\)。
第2節の proper marked triple-factor／rank 閉条件を適用すれば、
CA 全 rank四閉包の各点には \(\ell^3\mid D\)、\(\operatorname{rank}N_\ell\le1\) の label が残る。

一方 M1／M11 では \(\ell^3\mid D\) の唯一 hyperplane は h。
その homogeneous modular 場は全 source で
\(N_h=h\partial_x+x\partial_y\)、matrix rank 二。
従って両25像とも CA に含まれない。
CA が28次元であることだけを理由に包含を排除したのではない。

## 4. CB の全 magnetic source と二つの double label

rawCB の定義も import／実行せず読解し、全十六 cubic を保持した。
generic binary f は actual affine 変更で
\(f=X(h-2X)/(6a)\)、\(a\ne0\) と取れる。
\[
 P_{Xy}=f,\quad P_{zw}=\mu f,\quad
 P_{yz}=-f_Xz/2-H_{3,w},\quad P_{yw}=-f_Xw/2+H_{3,z},
 \quad P_{Xz}=P_{Xw}=0.
\]
全 divergence は \(2f_X\partial_y\)。
row cofactor \(P^{\ell j}/\ell\) は
\[
 Y_X=(h-2X)\partial_y/(6a),\qquad
 Y_{h-2X}=-2X\partial_y/(6a).
\]
どちらも trace零であり、factor-Casimir lift は \(P-E\wedge Y_\ell\)。
従って全 cubic・全 \(\mu\) に対して
\[
 N_X=2f_X\partial_y-5Y_X=(-3h+2X)\partial_y/(6a),
\]
\[
 N_{h-2X}=2f_X\partial_y-5Y_{h-2X}=(2h+2X)\partial_y/(6a).
\]
両 homogeneous matrix とも rank 一以下。
第2節の proper double-factor／rank incidence により、CB 全 rank四閉包にもこの label 条件が残る。

M1 の一般 \(D=h^4L\) では二乗 factor となる hyperplane は h ただ一つ。
従ってどちらか一つの double label を保持するだけでも \(N_h\) rank 一以下を強制し、
target の rank 二に矛盾する。
両 label を保つなら二つとも h へ coalesce し、simple label は L へ行く。
coalescence で rank 閉条件が失われることはない。

M11 の一般 \(D=h^3L_1L_2\) は、CB の proper factor incidence
\(D=H_0L^2M^2\) 自体の外。
二 label が同じ h なら四乗を要求し、異なるなら simple \(L_i\) の二乗を要求して、いずれも不可能。
従って M11 では追加の rank 条件も不要である。

## 5. 全旧・新25成分と二者の相違

rank二の全六成分は closed rank 条件で除外、rank四 dimension二十四以下は正確な25下界で除外する。
残る旧・新25次元成分への包含は、閉既約像の同次元性により一致を意味する。

M1 は全 \(D=H^4L\) の proper 閉条件内で、全 prior129、非定数g新成分、
新 constant-g 二共鳴と scalar-cubic の各 generic D はこの条件の外。
M11 は all-linear かつ triple-factor の閉条件内。
MPS25 の reduced D、positive NEW25 の rank五 Q、非定数gの rank五 Q、
constant-g 二共鳴の irreducible cubic はその generic 型と異なる。
SC は最大 multiplicity 二なので triple-factor の外。
nil31-DDH は最大 multiplicity 二、half／third の全 D は \(3H+2L\) の proper 閉条件内。
generic M11 の \(3H+L_1+L_2\)、\(L_1\ne L_2\) はその外である。
従って全既知25との一致は起こらない。

M1 対 M11 も、M1 の四乗 factor 閉条件に M11 の generic D が入らないので異なる。
七つの g零 B1–B6／V1 の新入口との相違も、両 M 像の all-linear 条件、
又は V1 の最大 multiplicity 二を用いた先行七族証明と一致する。
両 M source に交わる特別な閉境界が無い、と主張したものではない。

## 6. 確認した範囲と次の局所義務

元 source の全係数、全 \(\lambda\)／affine-v、全 F、actual stabilizer／scalar／係数逆、
unique repeated-H、CA／CB の全磁気と proper marked 閉条件、全 prior と相互識別を独立に確認した。
旧 rawCA／CB code は定義を読んだだけで一切実行・import・変更していない。

**成分性と未知の高次元親は未証明。**
両 whole source の指定した各一 fixed point で fullambient rank100 を得たとき、
この25次元下界と Euler24／scalar 一を除く接上界が一致するという最弱の矢印は正しい。
本監査ではその新しい点の選択又は演算を行っていない。
全 global catalogue／component count／主稿／Git は integrator の別責任として変更しない。
