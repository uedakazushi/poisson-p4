# 重複超平面層：低階数 modular 場と半単純階数二の全被覆

2026-11-07。主課題 E_multiple_hyperplane。
重複超平面を無限遠にした全二次 affine Poisson 構造
`P=C+L+Q`、`Pf(Q)=0` を扱う。原点の C は退化していてもよい。
modular 場は係数規約
`X^j=Σ_i ∂_i P^{ij}` により `X=c+Ax` と書く。
`tr A=0,[X,P]=0,Q∧Ax=0` を用いる。
全次数 Jacobi と `Pf(Q)=0` を保持し、任意の Q を Q45 と同一視しない。
第2節は二人が別々に全 Jacobi・座標 shear・閉包の境界を手で再導出し、
第1・3節は例外層担当者が独立に証明を読んで受理した。
第4–5節の普遍 Jacobi・modular・Pfaffian・全係数の座標恒等式は
`reports/CLASSIFICATION_RANK1_MODULAR_SOURCES_20261107.json` で QQ 上検算済み。
これは producer 自身の記号検算であり、別作者による独立再実装とはしない。
第4.1–4.2節の全CA/CB認識、任意三次項、b0=0境界への閉包は、
別担当が座標変換・Euler補正・全係数を手で再導出して独立に受理した。

## 1. rank A≤2 の有限な affine 場の型

rank A≤2、tr A=0 の Jordan 型は、次の五つに限られる。

| A | affine 平行移動・中心化群で残す c の代表 |
|---|---|
| 0 | 0、非零定ベクトル |
| 非零 nilpotent、型 (2,1,1) | 0、長さ2鎖の top、長さ1鎖の非零元 |
| nilpotent、型 (3,1) | 0、長さ3鎖の top、長さ1鎖の非零元 |
| nilpotent、型 (2,2) | 0、長さ2鎖の一つの top |
| diag(λ,−λ,0,0)、λ≠0 | 0、零固有空間の非零元 |

平行移動は c を `im A` の元だけ変えるので、まず `V/im A` の類を残す。
nilpotent の中心化群は、長い鎖の top から短い鎖の top への成分を加えられる。
非零類に現れる最長鎖の長さを選べば、同じ長さの GL 作用とこの加算で上表の代表になる。
異なる長さから最長鎖への加算はできないため、表の枝を混同しない。
半単純型の零固有空間には GL2 が作用する。
従って c を無断で零にすることなく、十二の affine 型を保持できる。

## 2. 階数二半単純型の全支持は MPS40 の和に入る

**命題。** `A=λ(x∂x−y∂y)`、λ≠0 の全一般 rank4 Poisson 構造は、
既知 MPS40 成分閉包の和に含まれる。c は `im A` に入る場合と入らない場合を全て含む。

### 2.1 全方程式を保った小さい表示

座標を `(x,y,z,w)` とし、平行移動の後で

\[
 X=\lambda(x\partial_x-y\partial_y)+c_0,
 \qquad c_0=0\ \hbox{または}\ \partial_z
\]

とする。二次以下の多項式二ベクトル上で、半単純作用と定数場の作用は可換である。
非零重み空間上では「非零スカラー＋nilpotent」が可逆なので、
`[X,P]=0` は半単純重み零と `[c0,P]=0` の両方を与える。

重み零の全係数を並べ、さらに `Q∧Ax=0` と modular の全係数を比較すると

\[
\begin{aligned}
 P^{xy}&=a xy+f(z,w),& P^{zw}&=g(z,w),\\
 (P^{xz},P^{xw})&=xU,&
 (P^{yz},P^{yw})&=y(c_0+X_g-U),                 \tag{1}\\
 X_g&=g_w\partial_z-g_z\partial_w,&
 a+\operatorname{div}U&=-\lambda.
\end{aligned}
\]

ここで f は次数高々二、g は affine-linear、U は affine-linear 場である。
特に `Q∧Ax=0` は、Pzw の二次項と xy 項を零にし、交差二次係数を反対符号にする。
定数項・一次項を落としていない。c0=∂z の場合は全係数が z に依らない。
以下で(1)に残る **全** Jacobi 条件を課す。

### 2.2 c0=0

`x≠0` で `t=xy` を使うと

\[
 P=x\partial_x\wedge Y+\operatorname{Jac}(tg),\qquad
 Y=(at+f)\partial_t+U,
\]

となる。三変数 `(t,z,w)` の Jacobian 二ベクトルは
`{t,z}=t g_w,{z,w}=g,{w,t}=t g_z` の意味である。
Jacobi は `[Y,Jac(tg)]=0`、従って

\[
 U(g)=(\operatorname{div}U)g,\qquad fg\in\mathbb C             \tag{2}
\]

と同値になる。これは多項式恒等式の比較なので、x=0 の点を捨てない。
g=0 では rank≤2 となるため除かれる。

**g が非定数の場合。** (2)から f=0。affine 座標で g=z とすれば

\[
 U=\alpha z\partial_z+(b_0+bz)\partial_w,\qquad a=-\lambda-\alpha.
\]

α≠0 の稠密開では `w↦w−(b/α)z` により b を消せる。
その後 P は、三変数 `(x,y,z)` の対数型二ベクトル Q3 と、
それを保つ diagonal な一次場 Y により `P=∂w∧Y+Q3` と書ける。
従って

\[
 P_\varepsilon=(1+\varepsilon w)\partial_w\wedge Y+Q_3
\]

は全 ε で二次 Poisson である。ε≠0 で w を平行移動すれば全て対数型となり、
Pfaffian は非零定数倍の xyzw、射影退化因子は hxyz w で SNC。
従ってこの pencil と P は MPS 成分に入る。
α=0 も、元の全パラメータ `(λ,α,b0,b)` の Poisson 族の閉包で保持される。

**g=g0≠0 が定数の場合。** (2)から f=f0 は定数、div U=0。
U の線形部分が可逆かつ半単純である稠密開では平行移動で定数部分を消し、
二変数線形変換で `U=μ(z∂z−w∂w)` とできる。
この場合、次の全パラメータ族は二次 Poisson である。

\[
 (-\lambda xy+f)\partial_x\wedge\partial_y
 +(x\partial_x-y\partial_y)\wedge\mu(z\partial_z-w\partial_w)
 +(g+\varepsilon zw)\partial_z\wedge\partial_w.                 \tag{3}
\]

f=g=0、ε≠0 なら Pfaffian は `−λεxyzw`、従って SNC 点となる。
元の f0,g0,ε=0 の点も、同じ既約な実 Poisson 族に含まれるので MPS 閉包に入る。
全 trace 零 affine 場 U のパラメータ空間で上の可逆半単純開は稠密である。
既知40閉包の有限和は閉なので、特異・nilpotent 線形部分と全定数部分も被覆される。

### 2.3 c0=∂z

全係数は z に依らない。(1)を

\[
\begin{gathered}
 g=g_0+d w,\quad U=A(w)\partial_z+B(w)\partial_w,
 \quad a+B'=-\lambda,\\
 P^{yz}=y(1+d-A(w)),\qquad P^{yw}=-yB(w)
\end{gathered}
\]

と書く。A,B は affine 一次式、f は w の次数高々二の式である。
四つの Jacobi 成分を展開すると、残る条件は正確に

\[
 B g'-gB'=0,\qquad g f'+(1+g')f=0.             \tag{4}
\]

**g≠0 の場合。** g が非定数なら `B=ρg`。
w の平行移動により g=dw、d≠0 とすると、第二式の全解は f=0、または

\[
 f=cw^n,\qquad d=-1/(n+1),\qquad n=0,1,2.     \tag{5}
\]

これらの共鳴枝を一つも捨てない。
Pfaffian は `g((ρ−λ)xy+f)` であり、これが非零となる開だけを使う。
ρ≠0 の稠密開では `z↦z−A_1w/(ρd)` によって A の w 係数を消せる。
その後

\[
 P=\partial_z\wedge Y+Q_3,
 \quad Y=-A_0x\partial_x+(A_0-1-d)y\partial_y+d w\partial_w,
\]

となり、Q3 は `(x,y,w)` の二次 Poisson 構造で `[Y,Q3]=0`。
この最後の条件は f 項についてまさに (5) であり、f=0 の場合は自動的である。
従って `Pε=(1+εz)∂z∧Y+Q3` は全 ε で二次 Poisson。
各共鳴族でも c=0 とでき、さらに ρ≠λ、ε≠0 とすれば
全係数が対数型で Pfaffian は非零定数倍の xyzw となる。
全 `(ρ,λ,c,A0)` の既約族とその閉包により、ρ=0 や ρ=λ の元も保持される。

g が非零定数なら (4) は B 定数、f=0を強制する。
これは、正規化前の `g=g0+dw,B=ρg,f=0` という実 Poisson 族で
d≠0 を稠密とする閉包に含まれる。w 平行移動の分母 d を境界まで延長する必要はない。

**g=0 の場合。** (4)から f=0、一般 rank4 条件は B≠0 である。
B の一次係数が非零の稠密開では、w を平行移動して `B=B1w` とし、
`z↦z−A1w/B1` により A を定数にできる。
再び P=∂z∧Y+Q3、Y は diagonal 一次場、Q3 は対数型である。
同じ `1+εz` pencil の一般元は SNC。
B1=0 の非零定数 B も、元の全 `(A,B,λ)` 族の閉包に含まれる。

各 pencil は実際の全 Poisson 族であり、SNC 点を含むその既約像閉包は
SNC locus と交わる全成分を尽くす MPS 定理を適用できる。
既知点への単なる形式退化だけで全支持を推論していない。これで命題が従う。□

## 3. A=0 についての部分還元

**c=0。** modular 恒等式と一般 rank4 から `Pf(P)` は非零定数となる。
この枝は全二次多項式 symplectic 構造、すなわち退化五次式 h⁵ の全支持である。
一次項が零とは限らず、beta零層または Q45 に無断で帰着しない。

**c≠0。** c=−κ∂t、κ≠0 と正規化する。
`[c,P]=0` より全係数は t に依らない。
三変数を x とすれば、全解は

\[
 P=\partial_t\wedge Y(x)+\operatorname{Jac}(f(x)),\qquad
 \operatorname{div}Y=\kappa,\quad Y(f)=\kappa f,
 \quad\deg Y\le2,\quad\deg f\le3                         \tag{6}
\]

へ還元される。三変数 bivector の divergence 零から一形式が閉となり、
多項式 Poincaré 補題で f を得る。f の加法定数は κ≠0 により一意に調整できる。

f の三次最高項 f3 が滑らかな平面三次曲線を定める場合は、(6)の全支持を被覆できる。
その三偏微分は正則列なので、二次 syzygy `Y2(f3)=0` は Koszul syzygy であり、
Y2 は `Jac(f3)` に関する一次関数の Hamiltonian 場である。
affine shear `t↦t+ℓ(x)` により、全 f の Hamiltonian 場を引いて Y2 を消せる。
残る一次部は `Y1(f3)=κf3` を満たす。
次数一の syzygy は存在しないため `Y1=(κ/3)Euler`。
平行移動で Y0 を消すと f は純三次斉次式になる。

この全正規形は

\[
 P_\varepsilon={\kappa\over3}(1+\varepsilon t)
                   \partial_t\wedge\operatorname{Euler}_3
                  +\operatorname{Jac}(f_3)
\]

という実二次 Poisson 族に入る。全三次式 f3 を自由に動かせ、
f3=xyz、ε≠0 で SNC 点を持つため、全正規形は MPS40 の和に含まれる。

f3 が特異、または deg f≤2 の全(6)では、追加 syzygy を残す必要がある。
この枝と、第一節の三 nilpotent Jordan 型の全支持被覆は未完了である。

## 4. rank一 nilpotent、c=0 の全源と二つの被覆済み枝

`X=x∂y` とする。`Q∧X=0` より `Q=∂y∧Z2`。
さらに modular の二次最高項比較から Z2 は y に依らず、div Z2=−x となる。
`[X,P]=0` の y を含まない成分は、残る affine 二ベクトルの y 依存を消す。
従って全点は、q=(x,z,w) を使って

\[
 P=\partial_y\wedge(Z(q)+y b)+R(q),\quad
 \deg Z\le2,\quad\deg R\le1,\quad b\in\mathbb C^3              \tag{7}
\]

と書ける。全 modular 条件と X 不変条件は

\[
 b=-\operatorname{div}R,\quad \operatorname{div}Z=-x,
 \quad R^{xj}=x b_j.                                           \tag{8}
\]

特に bx=0。B=bz∂z+bw∂w とすると、(8)の全解は

\[
 R=x\partial_x\wedge B+
 (\ell_0+\ell_x x+2b_z w-2b_w z)\partial_z\wedge\partial_w.     \tag{9}
\]

(9)は全パラメータで Poisson である。(7)の全 Jacobi は正確に

\[
 [R,R]=0,\qquad \mathcal L_bR=0,\qquad
 \mathcal L_ZR+Z\wedge b=0.                                   \tag{10}
\]

実際、y を含まない三添字が第一式、y を含む三添字の y 係数と定数係数が
第二・第三式を与える。従って低次 Jacobi を省略した最高項だけの表示ではない。

### 4.1 b=0、g定数：全源が CA に入る

このとき `R=g(x)∂z∧∂w`、g=ℓ0+ℓx x。
g=0 は rank≤2 なので、以下 g≠0 とする。
(10)と div Z=−x は

\[
 Z^x=F(x),\quad
 \operatorname{div}_{z,w}Z=F g'/g,\quad (Fg)'=-xg              \tag{11}
\]

に同値である。g=g0≠0 が定数なら

\[
 Z=(K-x^2/2)\partial_x+X_H,
 \qquad X_H=-H_w\partial_z+H_z\partial_w,\quad\deg H\le3.
\]

H には x のみの多項式を加えてよい。全 H を斉次化すると、(7)は

\[
 (x^2/2-Kh^2)\partial_x\wedge\partial_y+
 g_0h^2\partial_z\wedge\partial_w+\partial_y\wedge X_{H_3}
\]

となる。これは既存 CA の実源そのものである。非零 g0 は線形座標の尺度で吸収できる。
依存する普遍源は `sources/v06.zip` 内 `ambient_families.tex` の式 (PA)、
現在の全源実装は `code/ambient_CA_self_pair_all_nonzero_flags_20261103_v63.py` の rawCA。
その全空間成分性は `KNOWN_COMPONENT_AUDIT.md` の既存監査を用いる。

### 4.2 b=0、g非定数：全源が CB に入る

g=ax+b0、a≠0 とする。x=0 は modular 場の零超平面なので、g の根へ移して捨てない。
(11)を積分し、F が多項式となる条件を課すと

\[
 F=-\tfrac16(2x-b_0/a)(x+b_0/a),\quad
 \nu=Fg'/g=-x/3+b_0/(6a),\quad
 Z=F\partial_x+\tfrac\nu2(z\partial_z+w\partial_w)+X_H.       \tag{12}
\]

ここでも H は任意の次数高々三の多項式である。
斉次一次式 G=ax+b0h、L=2ax−b0h とすると斉次化した全係数は

\[
 P^{xy}=GL/(6a^2),\quad P^{zw}=hG,\quad
 P^{yz}=-Lz/(12a)-H_{3,w},\quad
 P^{yw}=-Lw/(12a)+H_{3,z},\quad P^{h\bullet}=0.                \tag{13}
\]

b0≠0 の開では H0=L、X0=G を新しい独立座標とする。
この座標変更の後、Euler 鉛直二ベクトル
`E∧((X0/(3a))∂y)` を引けば H0 が Casimir となる。
残る係数は正確に

\[
 \begin{gathered}
 f=P^{X_0y}=X_0(H_0-2X_0)/(6a),\qquad
 P^{zw}=-(2a/b_0)f,\\
 P^{yz}=-\tfrac12 f_{X_0}z-\widetilde H_{3,w},\qquad
 P^{yw}=-\tfrac12 f_{X_0}w+\widetilde H_{3,z}.                \tag{14}
 \end{gathered}
\]

従って、非零の比例定数を z,w の尺度で吸収した CB の実源である。
普遍源は上記 v06 の式 (PB)、比例係数を保持する現行実装は
`code/ambient_CA_self_pair_CB_source_rank_20261103_v63.py` の rawCB(f,mu*f,H)。
b0=0 の点も、元の正規形(12)の既約パラメータ源
`{a≠0}×A1_b0×(全H)` の閉包に含まれるので、全て CB 閉包に入る。
座標変換の分母 b0 を零へ代入する議論ではない。

### 4.3 b≠0：全10パラメータ源と新25次元成分への被覆

(z,w)の線形変換で b=∂z とする。(9)の二定数は
`w↦w+ℓ0/2+(ℓx/3)x` という affine shear で消せる。
よって

\[
 R=\partial_z\wedge V,\qquad V=-x\partial_x+2w\partial_w.
\]

(10)の ∂x∧∂w 成分は `2w∂zZx+x∂zZw=0` である。
次数条件により `Zx=A xz+F(x,w), Zw=−2A wz+G(x,w)`。
残る二成分は、`Zz=A z²/2+zD(x,w)+K0(x,w)` と書いて

\[
 (V+2)F=xD,\qquad (1-V)G=2wD,\qquad F_x+G_w+D=-x.          \tag{15}
\]

単項式 x^i w^j の V 重み −i+2j を比較すれば、全解は

\[
\begin{aligned}
 Z^x&=x(Az+d_0+f x+d w/3),\\
 Z^w&=w(-2Az-2d_0+(-1-2f)x-2d w/3),\\
 Z^z&=A z^2/2+z(d_0+d w)+K_0(x,w),\qquad\deg K_0\le2.       \tag{16}
\end{aligned}
\]

A,d0,f,d と K0 の六係数を合わせた A10 が、正規化した全支持を全射的に与える。
全パラメータで厳密に Poisson、modular 場 x∂y、Pfaffian は標準順序 (x,y,z,w) で

\[
 \operatorname{Pf}(P)=-(2wZ^x+xZ^w)=x^2w.
\]

従って一般階数四は全源で成立し、実射影因子は h²x²w。
この全源は、`RANK2_NIL_31_COMPONENT_REVIEW.md` で独立確認した新25次元成分に入る。
以下の橋は例外層担当者が提案し、本ノートの担当者が全 cofactor と支持条件を再導出した。

正準 lift `Π=P−E∧(x∂y)/5` に対し、もう一つの二重超平面 x=0 の cofactor は

\[
 Y_x=\Pi^\sharp dx/x
 =-(Az+d_0h+(f+1/5)x+d w/3)\partial_y+h\partial_z.
\]

従って `Π−E∧Yx` は x を Casimir にし、x=1 の chart における modular 場は

\[
 (5Az+5d_0h+5d w/3+5f+1)\partial_y-5h\partial_z.             \tag{17a}
\]

A≠0 の稠密開では線形部は rank二 nilpotent (3,1)、定数部はその像に入る。
平行移動と線形 Jordan 正規化後、その全15源に入る。
chain の末端一次式 `im(N*)²` は h の線であり、w は h と独立な不変一次式である。
新 affine Pfaffian は非零定数倍の h²w なので、この末端一次式だけの多項式ではない。

一方、新25成分の入口の全15源で n7=0 とすると、全 Jacobi が
`n11*n5=n11*n6=0` を強制する。その Pfaffian は

\[
 x_3(n_5+n_6x_3)(n_{10}x_3+n_9)
\]

となり、必ず末端一次式 x3 だけの多項式である。
従って(17a)から得た点は n7≠0 の全被覆済み開に入る。
全 A10 源は既約であり A≠0 は稠密なので、A=0 の全境界も新成分の閉性から含まれる。
これにより `X=x∂y,c=0` の全一般 rank四支持は、CA、CB、この新25成分で被覆された。
c が im A の外にある二枝は次節に残す。

## 5. rank一 nilpotent、cがim Aの外にある二枝の全還元

`X=x∂y+∂x` と `X=x∂y+∂z` をそれぞれ

\[
 (t,r,s,w)=(x,y-x^2/2,z,w),\qquad
 (t,r,s,w)=(z,y-xz,x,w)                                      \tag{17}
\]

で書き換えると X=∂t になる。この変換は Jacobian 一の多項式自己同型だが
**射影線形ではない**。従って、変換後の構造の既知成分への帰属から元の帰属を推論しない。
以下は全低次数方程式の正確な表示を得るためだけに使う。

両枝とも全点は、q=(s,w) 上の affine 場 U、affine 式 g、
次数高々二の F、次数高々三の H により

\[
 \begin{aligned}
 P'&=\partial_t\wedge Y+\operatorname{Jac}(f),\\
 Y&=-(F+(1+d)r)\partial_r+U,
 \quad d=\operatorname{div}U,\quad f=rg+H,\\
 U(g)&=d g,\qquad (U+1)H=Fg                                \tag{18}
 \end{aligned}
\]

と書ける。Jac は (r,s,w) 順の三変数 Jacobian 二ベクトルである。
H の加法定数は第二式に合うよう一意に決める。

全支持との同値を明示するため、逆変換後の六係数を記す。
第一の枝では U=u∂s+v∂w、r=y−x²/2、s=z として

\[
\begin{aligned}
 P^{xy}&=-F-(1+d)r,&P^{xz}&=u,&P^{xw}&=v,\\
 P^{yz}&=xu+r g_w+H_w,&P^{yw}&=xv-rg_s-H_s,&P^{zw}&=g.
\end{aligned}                                                   \tag{19}
\]

第二の枝では r=y−xz、s=x として

\[
\begin{aligned}
 P^{xy}&=-xu-rg_w-H_w,&P^{xz}&=-u,&P^{xw}&=g,\\
 P^{yz}&=F+(1+d)r-zu,&
 P^{yw}&=xv+zg-rg_s-H_s,&P^{zw}&=v.                            \tag{20}
\end{aligned}
\]

全係数は元の座標で次数高々二である。
逆に第4節冒頭と同じ最高項・modular・不変条件の比較を c≠0 で行うと、
それぞれ(19)、(20)の全係数を得る。変換後の modular 場が ∂t であることは
div Y=−1 と等価で、三変数部分は divergence 零なので Jac(f) と積分できる。
残る全 Jacobi は `Y(f)=−f`、その r 係数と定数係数が(18)である。
従って補助的な非線形座標で現れうる任意高次数項を切り捨ててはいない。

さらに g で次の有限な分岐に分かれる。

* g が非定数なら、q の affine 座標で g=s とできる。
  `U=αs∂s+(b0+bs)∂w` となり、第二式を s=0 に制限すると
  `(1+b0∂w)H(0,w)=0`。この作用素は有限次数の多項式上で可逆だから H=sK。
  従って `deg K≤2,F=(U+α+1)K` で全解を尽くす。
* g=g0≠0 が定数なら div U=0。
  H の最高三次項は `(U1+1)H3=0` を満たし、F=(U+1)H/g0。
  U1 が nilpotent なら H3=0。半単純固有値が ±μ なら、H3≠0 の可能性は
  `μ∈{±1,±1/3}` に限られる。これは全四つの三次重み ±μ,±3μ の比較である。
  特殊値も除外せず保持する。
* g=0 なら F は任意の二次式で、残る条件は `(U+1)H=0`。
  一般 rank四条件は f=H が零でないことである。

これらは全元方程式と同値な有限次数の源であり、まだ全源の既知成分への包含を示さない。
特に(17)を PGL 作用の代わりにして被覆済みと数えない。

### 5.1 固定した両modular場の源は12次元以下

g が非定数なら、その全 affine 三係数を戻した基底は三次元。
各 g の fiber は上の U の三係数と K の六係数で尽くされるため、合計12次元以下である。
g=s という正規化の九次元源だけを、元の全固定場源の次元と混同しない。

g=g0≠0 では trace零 affine 場 U は五次元、H の二次以下部分は六次元である。
一般 U の H3 は零だから g0 の一次元を加えて12。
H3≠0 の場合、U1 の固有値 ±μ は `μ∈{±1,±1/3}` に限られる。
U の次元は四へ落ち、H3 の固有空間は一次元だけなので、やはり `4+6+1+1=12`。
nilpotent U1 では U1+1 が三次斉次式上で可逆であり、H3=0 となる。

g=0 の場合は F の六係数が自由である。
残る incidence `I={(U,H): U(H)=−H, 0≠H, deg H≤3}` が六次元以下であることを示せばよい。
U の affine Jordan 型によって次の表で尽くす。

| U1 の型 | U の次元上界 | H の固有空間 | 合計上界 |
|---|---:|---:|---:|
| 相異なる非零固有値、一つの単項式共鳴 | 5 | 1 | 6 |
| 相異なる非零固有値、二つの共鳴 | 4 | 高々2 | 6 |
| 固有値 λ,0、定数部が im U1 内 | 3 | 高々3 | 6 |
| 固有値 λ,0、零固有方向の定数部が非零 | 4 | 1 | 5 |
| U1=−I/k、k=1,2,3 | 2 | k+1 | 高々6 |
| 非半単純な重複非零固有値 | 4 | 1 | 5 |
| 両固有値零 | — | 0 | 空 |

最初の二行では定数部を平行移動で消す。
各単項式 s^i w^j の条件は `iλ+jμ=−1`、1≤i+j≤3。
一つの条件は固有値二変数に一条件、異なる二つは独立な二条件となる。
λ,μ が相異なる非零値なら、この三角形内の共鳴単項式は高々二つである。
三つ以上が同一直線上に並ぶ場合は、λ=μ または λμ=0 となり、別行に入る。

零固有値がある場合、非零方向を s とする。
定数部が像内なら λ=−1/k で H=s^k p(w)、deg p≤3−k。
最大は k=1 の三次元である。像の外の定数部が非零なら、その平行移動作用により
p'(w)=0 が必要となるので、固有空間は一次元だけとなる。
scalar の場合は中心を平行移動した後の斉次 k 次式全体が固有空間で、中心の二次元を戻す。
非半単純な重複固有値では λ=−1/k が必要であり、nilpotent 作用の kernel は
一つの一次式の k 乗だけである。最後の行の affine 場は多項式上で局所冪零なので、
固有値 −1 の非零多項式を持たない。

よって `dim I≤6`、F の六係数を加えて g=0 の全源も12次元以下となる。
全 Jordan・共鳴・零係数枝を含み、一般 U の固有空間次元を全点に適用してはいない。

### 5.2 全短鎖・長鎖の g=0 と定数gのβ像

`X=x∂y+∂z` の Aff4 安定化群は11次元である。
線形部の GL4 中心化群は10次元。
長さ一の鎖にある c の `V/im A` の類は二次元の軌道を持つため、
その固定は二条件であり、最後に ker A の三次元平行移動を加えると11となる。
従ってこの affine 場の軌道は `20−11=9` 次元である。
全固定場源の上界12を合わせ、全 affine 支持は高々21次元。
非退化な基点値を固定する評価 fiber では GL4 の六次元推移軌道を除き、
primitive β 像の次元は高々15となる。全固定 symplectic 値 J に同じ上界が成り立つ。
既存の perspective の橋により、元の重複超平面層でもこの型は W=A16 を優越しない。

一方 `X=x∂y+∂x` は、c の商空間軌道が三次元なので Aff4 軌道は10次元である。
単なる源次元の比較は `10+12−6=16` までしか与えないが、g=0 と定数 g は以下で閉じる。

**g=0 の全支持。** (19)で F は Pxy だけに現れ、全 Jacobi と Pfaffian から独立な任意二次式である。
任意の symplectic 基点の q 座標を q0 とすると、q−q0 の斉次二次式三次元を F に加える実族は
元の座標の全0/1 jet を固定する。従って全 affine 支持の primitive β 像は
`12+10−6−3=13` 次元以下。短鎖では同じ引算により12次元以下である。
ここでは実 affine fiber を使っており、一点の微分核から一般 fiber 次元を推定していない。

**g=g0≠0、H3=0 の全支持。** 全固定場源は
`Gm_g × A5_U × A6_H` であり既約12次元。F=(U+1)H/g で一意に決まる。
その Aff4 飽和は高々22次元、無限遠超平面の四次元を戻して射影化した
既約閉包 Z は高々25次元である。
この Z は旧 positive-weight の25次元成分 NEW25 に等しい。

包含を確かめる入口は `reports/Q_POSITIVE_WEIGHT6_NEW_FAMILY_20261103_v63/certificate.json.gz`
の `full_affine_P` である。座標を (x,y,z,w)、源の係数を n0,n1,n2 とすると、その modular 場は
\[
 Z_0=\frac{4n_0}{3}\partial_y+
       (n_0x+n_0n_2y+n_1/48)\partial_z.
\]
稠密開 n0 n2≠0 で d=4n0²n2/3 と置き、
\[
 X=(n_0x+n_0n_2y+n_1/48)/d,\qquad Y=z/d,\qquad q=(x,w)
\]
と線形変更すれば `Z0(X)=1,Z0(Y)=X,Z0(q)=0`。
従ってこれは長鎖の正規 modular 場で、`g={x,w}=8/9` は非零定数である。
元の Pfaffian は次数高々二であり、(18)の `f=g(Y−X²/2)+H(q)` と一致するので H3=0。
旧 NEW25 の全三変数源の稠密開とその PGL 飽和が Z に入る。
旧成分性・次元25を `NEW25_COMPONENT_AUDIT.md` の依存から入力すると、
既約閉集合 `NEW25⊂Z, dim Z≤25` より Z=NEW25。
したがって n0 n2=0 の旧境界だけでなく、現在の全 Gm×A11 源も同成分に入る。
この証明は一代表の所属を全源へ拡張していない。
旧成分の primitive 像が proper である既存結果も接続する。

**g=g0≠0 の三次共鳴源。** U1 の三次固有多項式 T≠0、`U1(T)=−T` を一つ固定する。
任意の symplectic 基点 q0 において
\[
 \delta H=T(q-q_0),\qquad \delta F=(U+1)\delta H/g_0
\]
と置く。δH は三次以上に消える。δF の最高三次項は `(U1+1)T=0` により消えるので、
δF は許容される二次式であり、q0 では二次以上に消える。
U と g0 を固定して H,F にこの方向を加えると、全 Jacobi を保つ実 affine line となる。
(19)の δPyz,δPyw は δH の一階微分、δPxy は −δF だから、全0/1 jet を固定する。
この実一次元 fiber と源次元上界12により、全 primitive β 像は `12+10−6−1=15` 以下。
H3 の係数が零の境界点でも同じ affine line が存在し、境界を除外しない。

§5.1 の全 Jordan 次元表と短鎖の15次元上界は generic-structure 担当が独立に手で確認した。
本節の g=0 の実三平面、定数 g の旧25次元成分への包含、三次共鳴の実直線は同担当の提案を
本ノート担当が全 jet と源次元の論理まで独立に確認した。
残る長鎖 g 非定数の全源は、次節の別の実有理曲線によって非優越となる。
自由 K の二次項は一般 symplectic 点で一次 jet に効くため、g=0 の三次元 fiber を流用しない。

### 5.3 長鎖・非定数gの全点を通る同じ一次jetの実曲線

**命題。** `X=∂x+x∂y`、g 非定数の全 rank四源の primitive β 像は次元高々15である。

まず q の affine 座標を使い g=z とする。このとき

\[
 U=az\partial_z+(b+cz)\partial_w,\qquad
 H=zK,\quad F=(U+a+1)K,\quad r=y-x^2/2,
\]

で K は任意の二次以下の多項式である。
元の六括弧は (19) に正確に戻り、Pfaffian は `−z(r+K)`。
従って任意の symplectic 基点では z≠0。
X を保つ変更 `x↦x−x0,y↦y−x0 x+constant`、z の非零 scalar 倍、w の平行移動を使えば、
基点を `(x,y,z,w)=(0,0,1,0)` に取れる。z の scalar 倍は g=z の表示も保つ。
さらに基点で零の affine 式を y に加えることで K の二つの一次 Taylor 係数を消し、

\[
 \zeta=z-1,\qquad K=k+A\zeta^2+B\zeta w+Cw^2,\qquad k\ne0
 \tag{20a}
\]

とできる。ここで行った変更は全て affine であり、非線形な r 座標を射影同値とみなしていない。

`v=b+c`、`ℓ=vζ−aw` と置く。λ≠0 に対して

\[
 a_\lambda=\lambda a,\qquad b_\lambda=b,\qquad
 c_\lambda=\lambda v-b,\qquad
 K_\lambda=\lambda K+\frac{\lambda(\lambda-1)}2\ell^2
 \tag{20b}
\]

とする。これらを同じ全源 (19) に代入した二次 Poisson テンソルを Pλ とし、
基点を固定する可逆 affine frame

\[
 g_\lambda(x,y,\zeta,w)
   =\bigl(x+(1-\lambda)\ell,\ y/\lambda,\ \zeta,\ w\bigr)
 \tag{20c}
\]

で `Qλ=(gλ)*Pλ` と置く。Q1 は元の全テンソルである。
全有理恒等式として

\[
 Q_\lambda(p)=Q_1(p),\qquad
 dQ_\lambda|_p=dQ_1|_p,
 \qquad
 \left.\partial_\lambda\operatorname{coeff}_{x^2}
       (Q_\lambda^{yw})\right|_{\lambda=1}=-\tfrac12
 \tag{20d}
\]

が成立する。従ってこの曲線は全ての a,b,c,k,A,B,C において非定数で、完全な0/1jetを保つ。
元の座標で全係数は次数高々二のままであり、全 Jacobi を満たす。
単に形式変形の有限次数が解けるという主張ではない。

恒等式の導出も短く記す。
λ=e^(−t) と形式的に微分すると、源の変化は
`δa=−a,δv=−v,δb=0,δK=−K−ℓ²/2`。
全0/1jetに対して、その変化は線形場 `ℓ∂x+y∂y` の pushforward 微分で補償される。
この線形場の係数は t に沿って明示的に変化し、積分後の有限行列は (20c) になる。
存在の根拠は微分方程式の形式解ではなく、(20b)–(20d) の直接の有限有理恒等式である。

§5.1 により固定場の全源は12次元、affine modular orbit は10次元以下なので、
全 affine 飽和源は22次元以下。固定した非退化定数値 J の開は六次元を引いて16次元以下である。
各点の元の source 座標へ戻って (20a)–(20d) を適用すれば、全 fiber に非定数曲線が存在する。
元の座標に戻す変換は曲線に依らないため、全0/1jetの固定と非定数性を保つ。
同じ定数値と一次項を持つ全曲線には、同じ Darboux frame と同じ有限 perspective を適用できる。
その可逆性により primitive β の fiber も一次元以上であり、像次元は `12+10−6−1=15` 以下。
これで rank一 nilpotent の全 affine 型について、一般 β を優越しないことが従う。
各特殊源の既知成分帰属をここで全て決定したとはしない。

公開恒等式証明書は `reports/CLASSIFICATION_RANK1_LONG_ACTUAL_FIBER_20261107.json`。
`code/classification_rank1_firstjet.py --curve` は全六括弧、全四 Jacobi、全30 zero/first jet残差、
次数上界、非定数係数 −1/2 を QQ(a,b,c,k,A,B,C,λ) 上で検算する。
有限契約は `jobs/rank1_long_actual_fiber.json`、実測1.062秒。
標本点の微分階数15は発見のための診断として別 report に残し、本命題の上界には使わない。
独立検算 `CLASSIFICATION_RANK1_LONG_ACTUAL_FIBER_REVIEW_20261107.json` は
別の Fraction/Laurent 疎実装で source を K,U から再構成し、affine pushforward の全六括弧、
全四 Jacobi、全30jet差、次数上界、非定数係数 −1/2 を確認した（0.047秒）。
別担当二名が全基点の affine 正規化、全実 fiber と次元の引算を独立に読解し、
統合担当も全 rank一 G 非優越を受理した。

## 6. A=c=0 の h五乗層：既存rank二分類の適用範囲

旧 `CLAIMS.json` の NORMAL_CLASSIFICATION は「5Hは対象外」、
Q_FIVE_COEFFICIENT_FAMILY は「全5H locusの分類ではない」と明記している。
従って旧証明の全 h五乗支持被覆を接続できたとはしない。

一般 rank四かつ divergence 零では、modular 恒等式から Pf(P)=k≠0。
四次元体積形式 vol により

\[
 \alpha=\iota_P\mathrm{vol}=J+\beta+\gamma
\]

を係数次数 0,1,2 に分けると、全問題は厳密に

\[
 d\beta=d\gamma=0,\quad J^2\ne0,\quad J\wedge\beta=0,
 \quad 2J\wedge\gamma+\beta^2=0,\quad
 \beta\wedge\gamma=0,\quad\gamma^2=0                         \tag{21}
\]

となる。逆に(21)なら α は全 affine 空間で非退化な閉二形式となる。
Pfaffian が定数なので、その対応する P は定数倍の逆 symplectic 形式であり Poisson。
この双方向の同値は、未知の γ を primitive な Q45 条件へ置き換えない。

γ が零でないとき、対応する斉次二次 Q は厳密 Poisson、divergence 零、rank二である。
新しい変数 h を加えた h-Casimir lift は既存 v14 の rank二分類に適用できる。
ただし `RANK_TWO_SUPPORT.md` の結論は複数の閉包への包含であり、
それだけから全 γ に特定の Jacobian 正規形を付けることはできない。
以下は、そのうち **因子的零点がなく安定な接層を持つ分岐**に限る実表示である。

この分岐の v14 の実表示は

\[
 dh\wedge\gamma=dL_1\wedge dL_2\wedge dG_3
 \quad\hbox{または}\quad
 dh\wedge\gamma=dL\wedge dQ_1\wedge dQ_2.
\]

第一の場合は L1,L2 の少なくとも一つが h と独立な一次 Casimir である。
第二の場合も L が h と独立なら同様である。
その一次式を h=0 へ制限して ℓ とすれば `dℓ∧γ=0`。
一次 Koszul 恒等式と相対多項式 Poincaré 補題により
`γ=dℓ∧df3`、f3 は斉次三次式となる。
第二の場合で L が h に比例すれば、h=0 への制限から直接
`γ=dq1∧dq2`、q1,q2 は斉次二次式を得る。
従ってこの安定分岐では、実際の全源を

\[
 \gamma=d\ell\wedge df_3\quad\hbox{または}\quad
 \gamma=dq_1\wedge dq_2                                     \tag{22}
\]

へ縮約でき、残る条件は(21)の全 β 方程式である。
非安定、因子的零点を持つ分岐、または各閉包の未指定の特殊点を(22)で覆ったとはしない。
v14 の安定証明は `RANK2_STABLE_REVIEW.md` で独立に読解済みの結果を入力とする。

## 有限恒等式の再現

```sh
python3 scripts/classification.py gate research/classification/jobs/rank1_modular_sources.json
.venv/bin/python code/classification_rank1_modular_sources.py --output /tmp/rank1-modular-sources.json
```

job は E_multiple_hyperplane、60秒・1536MiB。実測0.343秒。
全10源の Jacobi・divergence・Pfaffian、二つの非線形補助座標からの全六係数の逆像と
degree≤2、定数 modular 表示の全 Jacobi が(18)の二残差の Jacobian となる恒等式、
全20三次係数を保持した(13)→(14)の全十係数一致を確認した。
CA/CB の既存全成分証明は再計算せず、認識された実源への包含だけを新規に検算した。
