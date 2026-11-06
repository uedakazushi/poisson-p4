# 五つの相異なる超平面を実退化因子とする全支持

2026-11-07。E_reduced_nonnormal。
以下の正常形、Jacobi 全条件、逆二形式の符号、実変形を独立に手計算した。
有限計算や標本からの一般性の推測は使わない。

**結論。** 一般階数4の Poisson 構造の実退化因子が五つの相異なる超平面の和なら、
構造は T01 の全閉包に含まれる。
一次式が5次元を張る場合は既存の toric 認識を使う。
4次元を張る場合は第1–4節、3次元の場合は第5–7節の実変形で全支持を覆う。
2次元以下の場合は第8節で一般階数4が不可能と示す。
4次元の場合、
唯一の線形関係の支持が4本または5本の場合は、そのような構造は存在しない。
支持が3本の場合は以下の全正常形と実変形で覆う。

## 1. 関係の支持が4本以上なら rank4 は不可能

五つの一次式の関係空間は一次元である。
相異なる超平面なので、その唯一の関係の支持の大きさ m は3,4,5のいずれか。
関係に現れる二因子を h と l と選ぶと、残る三因子 x,y,z と h は独立である。
残る第五座標を w とし、h=1 で

\[
 l=a_0+a_xx+a_yy+a_zz,\qquad a_0\ne0
\]

と書く。非零な a_x,a_y,a_z の数は m−2 である。
既存の h-Casimir lift により、全 Poisson 構造は affine 四変数の次数2以下の
strict Poisson bivector P で表せる。各退化因子は Hamiltonian 場で保存されるため、
各 x,y,z,l の主イデアルは Poisson ideal である。

x,y,z の保存と次数条件から

\[
 \{x,y\}=Axy,\qquad\{x,z\}=Bxz,\qquad\{y,z\}=Cyz
                                                        \tag{1}
\]

となる。`{l,x},{l,y},{l,z}` は斉次二次式であり、同時に l で割り切れる。
定数項非零の affine 一次式 l で割り切れる斉次二次式は零である。
実際、商の定数項と一次項を順に比較すると商が零になる。
従って l に二つ以上の x,y,z が現れれば (1) の A,B,C は全て零である。
このとき P は `V∧∂w` の形となり、一般階数は高々2。m=4,5 は排除される。

## 2. 三因子の関係の場合の全正常形

m=3 なら因子と座標を尺度変更して `l=x−1` とできる。
同じ保存条件から A=B=0 であり、次数条件によって

\[
 \begin{gathered}
 \{x,w\}=\kappa x(x-1),\qquad\{y,z\}=Cyz,\\
 \{y,w\}=yf(x,y,z,w),\qquad\{z,w\}=zg(x,y,z,w),
 \end{gathered}                                             \tag{2}
\]

f,g は affine 一次式となる。一般階数4は `κC≠0` と同値である。
四つの Jacobi 式を全て使うと、x,y,w と x,z,w から
`f_w=g_w=0`、y,z,w から `y f_y+z g_z=0` を得る。
x,y,z の式は恒等零である。従って正確に

\[
 f=p_0+p_1x+p_2z,\qquad g=q_0+q_1x+q_2y
\]

となり、他のパラメータ条件はない。
affine shear

\[
 w'=w+(q_2/C)y-(p_2/C)z
\]

は x,y,z,l を変えず、両方の混合係数 p2,q2 を消す。
以下 w' を w と書く。この結果、全支持は

\[
 \boxed{(\{x,y\},\{x,z\},\{x,w\},\{y,z\},\{y,w\},\{z,w\})
 =(0,0,\kappa x(x-1),Cyz,y(p_0+p_1x),z(q_0+q_1x))}
                                                            \tag{3}
\]

で尽くされ、κ,C は非零、p0,p1,q0,q1 は任意である。

## 3. 逆二形式と全パラメータ変形

逆二形式の行列は `Ω=−P^{-1}` とする cofactor 規約を用いる。
(3)の逆形式は

\[
 \omega_0={1\over C}{dy\over y}\wedge{dz\over z}
 +{dx\over\kappa x(x-1)}\wedge
 \left(dw-{p_0+p_1x\over C}{dz\over z}
             +{q_0+q_1x\over C}{dy\over y}\right).             \tag{4}
\]

符号は4×4 cofactor によって (3) の六括弧全てへ戻る。
`l_t=x−1+tw` と置き、`X=dlog x,L_t=dlog l_t,Y=dlog y,Z=dlog z` と書く。
次の閉二形式を考える。

\[
 \begin{split}
 \omega_t={}&{1\over C}Y\wedge Z+{1\over\kappa t}X\wedge L_t\\
 &-{1\over\kappa C}\big((p_0+p_1)L_t-p_0X\big)\wedge Z\\
 &+{1\over\kappa C}\big((q_0+q_1)L_t-q_0X\big)\wedge Y.
 \end{split}                                                \tag{5}
\]

すべて定係数の対数一形式の外積なので `dω_t=0`。
見かけの t 分母は

\[
 {1\over t}X\wedge L_t={dx\wedge dw\over x l_t}
\]

で消え、(5)は t=0 で厳密に (4)へ収束する。
順序 `(X,L_t,Y,Z)` での定数交代行列 B_t の Pfaffian は

\[
 \operatorname{Pf}(B_t)={1\over\kappa Ct}
       +{\Delta\over\kappa^2C^2},\qquad
 \Delta=p_0q_1-p_1q_0.                                    \tag{6}
\]

従って全パラメータについて一般の非零 t で非退化である。
Δ=0 の場合も (6) の最初の項は消えない。

さらに逆行列を明示すれば、極限が二次 Poisson 空間内で存在することも直接分かる。
`d_t=κC+tΔ` と書くと、(5)の逆 bivector の六括弧は

\[
 \boxed{\quad {\kappa C\over d_t}
 \left(
 \begin{array}{c}
 t(p_0+p_1)xy\\
 t(q_0+q_1)xz\\
 \kappa x l_t\\
 Cyz\\
 y\{p_1x+p_0(1-tw)\}\\
 z\{q_1x+q_0(1-tw)\}
 \end{array}\right) \quad}                                \tag{7}
\]

である。全係数が次数2以下、`d_0=κC≠0` なので t=0 を含む代数的な開上で正則であり、
t=0 の値は (3)に等しい。閉で非退化な逆形式から一般 t の Jacobi が従い、
多項式恒等式として t=0 にも成立する。新しい余分なパラメータ条件を課していない。

## 4. 全射影構造と T01 への帰属

(7)を h-Casimir として次数2に斉次化する。
実退化五次式は

\[
 {\kappa^2C^2\over d_t}\,hxyz(x-h+tw).                    \tag{8}
\]

t≠0 では五一次式 `h,x,y,z,x−h+tw` は独立である。
また (5)は対応する四つの affine 比の対数微分に関する定係数非退化二形式なので、
その逆はこれらの座標で logcanonical な Poisson 構造である。
従って一般 t の全構造は T01 の toric 開部分に入る。
T01 は閉じた全空間成分であり、(7)は t=0 まで正則な実変形なので、(3)も T01 に入る。

この証明は一般の p,q だけを扱ってから特殊値へ推測したものではない。
任意の p0,p1,q0,q1 と非零 κ,C に対して (5)–(8) が成立するため、
関係の支持が3本の場合の全支持を覆う。
射影変換と shear を戻せば、元の階数4配置の全点について結論が成立する。

依存は、超平面因子に対する h-Casimir lift、各退化因子の Poisson 保存、
および五独立超平面の toric 開部分が T01 に入る既存認識だけである。
Q45 の分類や全 S61 の一般ファイバーの尽くしには依存しない。

## 5. 一次式が3次元を張る場合の線配置の制約

一因子 h と独立な二因子 x,y を選び、h=1 とする。
残る二因子と合わせた四つの affine 直線を `ℓ_i(x,y)` と書き、自由座標を z,w とする。
`{x,y}=Axy` であり、第三の相異なる直線の保存から A=0 が従う。
従って P の一般階数4は、二つの混合ベクトル場

\[
 X_z=\{x,z\}\partial_x+\{y,z\}\partial_y,\qquad
 X_w=\{x,w\}\partial_x+\{y,w\}\partial_y
\]

の行列式が非零であることと同値であり、実退化因子から

\[
 \det(X_z,X_w)=c\prod_{i=1}^4\ell_i,\qquad c\ne0.             \tag{9}
\]

符号は c に吸収する。右辺は z,w に依らない。
z=w=0 に制限した二場も (9) を満たす。
その係数を h,x,y に関して次数2に斉次化し、h 成分を零とした場を A1,A2 とする。
Euler 場 E と並べると

\[
 \det(E,A_1,A_2)=cF,\qquad
 F=h\prod\widetilde\ell_i.                                \tag{10}
\]

三場は全て F の各因子に接する。
以下の自由基底の主張は Cramer 則だけで証明できる。
任意の logarithmic 多項式場を一列に代入した行列式は、各直線の一般点で
三列が二次元の接平面に入るため、その直線式で割り切れる。
F は平方自由なので各 replacement determinant は F で割り切れ、
任意の logarithmic 場は E,A1,A2 の多項式係数線形結合となる。
従ってこれらは次数 `(1,2,2)` の自由基底である。

特に h 成分が零の次数1以下の logarithmic 場は零である。
次数1の場合、Cramer 則で A1,A2 の係数の分子は次数4、分母 F は次数5だから零。
残る E の定数倍も h 成分が零なので消える。
次数0の場合も同様である。
元の X_z,X_w の z,w 依存係数は、全次数2以下という条件から、
この次数1以下の affine logarithmic 場である。従って全て零となる。
これは自由座標依存を仮定で捨てるのではなく、(9)–(10)から排除する議論である。

E を使って A1,A2 から F の対数微分を差し引くと、二つの自由な次数2の
gradient syzygy を得る。従って `S=C[h,x,y]` における Jacobian 商の Hilbert series は

\[
 {1-3t^4+2t^6\over(1-t)^3}
\]

であり、射影 Jacobian scheme の長さは12である。
五直線配置では各 m 重点の Tjurina 数は `(m−1)²`。
従って

\[
 \sum_p(m_p-1)^2=12,\qquad
 \sum_p\binom{m_p}{2}=10.
\]

差を取ると、三重点はちょうど二つで、四重点以上は存在しない。
二つの三重点を作る三直線ずつは、五本しかないので少なくとも一線を共有する。
二線を共有すれば二つの交点が同じになるため、共有線はちょうど一つである。
その線を h=0 に選び直すと、残る四直線は二組の平行線になる。
affine 座標変更で正確に

\[
 F=h\,x(x-h)y(y-h)                                      \tag{11}
\]

の配置へ帰着できる。

## 6. 格子配置での全 Poisson 正常形

h=1 とし `f=x(x−1),g=y(y−1)` と書く。
各直線の保存と次数2以下の条件により、混合括弧行列は

\[
 \begin{pmatrix}\{x,z\}&\{x,w\}\\\{y,z\}&\{y,w\}\end{pmatrix}
 =\begin{pmatrix}a_1f&a_2f\\a_3g&a_4g\end{pmatrix},
 \qquad a_1a_4-a_2a_3\ne0.
\]

係数 a_i は自由座標も含めた定数でなければならない。
fiber 座標 (z,w) の定数 GL2 変換でこの行列を `diag(f,g)` にできる。
残る括弧 `{z,w}=B` について、x,z,w と y,z,w の Jacobi は
`f B_z=g B_w=0` を与える。従って B は x,y の二次以下の多項式である。
shear `z'=z+αy,w'=w+βx` は B を `B+αg−βf` へ変えるので、
二乗項を除いて全正常形を

\[
 P=f\partial_x\wedge\partial_z
    +g\partial_y\wedge\partial_w
    +(a+bx+cy+dxy)\partial_z\wedge\partial_w               \tag{12}
\]

とできる。a,b,c,d は任意であり、残る Jacobi は恒等零。
この際にも未検査の自由座標依存や特殊値を捨てていない。

## 7. 格子配置の全パラメータ toric 変形

(12)の逆形式は

\[
 \omega_0={dx\wedge dz\over f}+{dy\wedge dw\over g}
       -{(a+bx+cy+dxy)\,dx\wedge dy\over fg}.             \tag{13}
\]

第3節と同じ cofactor 規約なので、最後の符号は負である。
`L=x−1+tz,M=y−1+sw`、
`X=dlog x,U=dlog L,Y=dlog y,V=dlog M` とし、

\[
 \begin{split}
 \omega_{t,s}={}&{1\over t}X\wedge U+{1\over s}Y\wedge V\\
 &-a(U-X)\wedge(V-Y)-bU\wedge(V-Y)\\
 &-c(U-X)\wedge V-dU\wedge V.                          \tag{14}
 \end{split}
\]

これは閉形式であり、見かけの t,s 分母はそれぞれ
`X∧U/t=dx∧dz/(xL)`、`Y∧V/s=dy∧dw/(yM)` で消える。
t=s=0 で (13)に一致する。
定数期間行列の Pfaffian は

\[
 {1\over ts}+bc-ad.
\]

逆 bivector も全六括弧を明示できる。
`p=a+b+c+d,q=a+b,r=a+c,δ=bc−ad,D=1+δts` と置くと、
順序 `(xy,xz,xw,yz,yw,zw)` で

\[
 \boxed{\quad {1\over D}
 \left(
 \begin{array}{c}
 ptsxy\\
 xL\\
 tx(qM-py)\\
 sy(px-rL)\\
 yM\\
 aLM-rLy-qxM+pxy
 \end{array}\right)\quad}                              \tag{15}
\]

となる。全係数は二次以下、D(0,0)=1 なので原点で正則であり、
t=s=0 の最後の成分は正確に `a+bx+cy+dxy` となる。
Pfaffian と斉次化後の実退化五次式はそれぞれ

\[
 -{xyLM\over D},\qquad
 -{hxy(x-h+tz)(y-h+sw)\over D}.                         \tag{16}
\]

t=s=ε とすれば、一般の ε≠0 で五一次式は独立、期間行列も非退化である。
従って一般構造は T01 の toric 開部分に入り、正則な ε=0 極限は任意の (12)となる。
これにより一次式が3次元を張る全支持も T01 の閉包に入る。

## 8. 一次式が2次元以下の場合

2次元の場合は、一因子を h=0 とすると、残る四因子は
`x−a_i`、相異なる四つの定数 a_i で書ける。
各因子の保存により、任意の座標 j に対する `{x,j}` はこの四因子全てで割り切れる。
次数は2以下なので `{x,j}=0`。
x は Casimir となり、残る三次元の交代行列の階数は高々2である。
従って一般階数4の構造は存在しない。
1次元では相異なる五超平面という前提自体が成立しない。

以上で、五つの相異なる超平面を実退化因子とする全ての一般階数4構造が覆われる。
結論は T01 への支持包含であり、同じ点を通る他成分の不存在や精密な局所スキーム構造を
主張するものではない。
