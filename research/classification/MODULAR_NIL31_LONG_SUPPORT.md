# 冪零(3,1)型・長い鎖の定数項：全六源と二つの22次元成分

2026-11-07。主課題 E_multiple_hyperplane。
modular 場を `X=x2∂x1+x3∂x2+∂x3` とする全 affine 二次 Poisson 構造を扱う。
入口は `CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json` case1 の
全11源と全 Jacobi である。全90線形係数還元は当該証明書を入力とし、ここで再実行しない。
本ノートはその全残条件を解き、二源を既存成分へ認識し、A源を実toric変形で被覆し、
さらに二つの新全成分を得る。最後の r=−1 の B源も、§5の全パラメータ実変形で
既済の三相異超平面＋既約二次因子層へ被覆する。これで全六源の帰属が閉じる。
第1–2節の全六源・全rank四条件・超平面変更による二源の全認識は、別担当が手で独立確認した。

## 1. 全六源の導出

入口の11変数を `n0,...,n10`、r=n10 とする。
n2,n3,n7,n9 は全ての枝で自由であり、全 Jacobi の残条件は

\[
\begin{gathered}
 n_5(r+1)=0,\quad n_1(3r+1)=0,\quad
 2n_1n_9+(2r+1)n_6=0,\\
 (r+1)n_0-n_3n_8+n_5n_7+n_6n_9=0,\quad
 n_8(3r-2)=0,\quad (2r-3)n_4-5n_5n_9=0.                 \tag{1}
\end{gathered}

r が `−1,−1/3,−1/2,2/3,3/2` の外なら、(1)は順に
n5=n1=n8=n6=n4=n0=0 を強制し、実 Pfaffian が零となる。
従って rank四ではこの五値だけを残せばよい。
相異なる五値を代入すると、除算で特別な枝を失うことなく次の全六源を得る。

| 名称 | r | 五つの自由変数 | 他の変数への全代入 | rank四条件 |
|---|---|---|---|---|
| A | −1 | n0,n2,n3,n7,n9 | n1=n4=n5=n6=n8=0 | n0≠0 |
| B | −1 | n0,n2,n3,n5,n9 | n1=n6=n7=n8=0,n4=−n5n9 | (n0,n5)≠(0,0) |
| C | −1/3 | n1,n2,n3,n7,n9 | n4=n5=n8=0,n6=−6n1n9,n0=9n1n9² | n1≠0 |
| D | −1/2 | n2,n3,n6,n7,n9 | n1=n4=n5=n8=0,n0=−2n6n9 | n6≠0 |
| E | 2/3 | n2,n3,n7,n8,n9 | n1=n4=n5=n6=0,n0=3n3n8/5 | n8≠0 |
| F | 3/2 | n2,n3,n4,n7,n9 | n0=n1=n5=n6=n8=0 | n4≠0 |

r=−1 では残る積条件が n5*n7=0 であり、A と B はその両枝である。
交わりを除外せず、片方を他方の稠密部分とみなしていない。
各源の指定された rank四開は既約である。

座標を `(x,y,z,w)=(x1,x2,x3,x4)` と略記すると、実 affine Pfaffian は順に

\[
\begin{array}{c|l}
A&n_0(n_9-w)\\
B&(n_9-w)\{n_0+n_5(-y+z^2/2)\}\\
C&-n_1(w-3n_9)^3/3\\
D&-n_6(w-2n_9)^2/2\\
E&-n_8\{x-yz+z^3/3+n_2+(3n_3/5)(w-n_9)\}\\
F&n_4(y-z^2/2+n_7).
\end{array}                                                   \tag{2}
\]

これにより表の rank四条件は必要十分である。
下記の四源 A,B,E,F は、四 raw Jacobi と実 Pfaffian を有限恒等式として検算した。
C,D は(1)への代入と次節の全支持帰着で扱う。

## 2. CとDは全て既存の冪零Jordan4成分に入る

五変数への斉次化後、全 modular 場は

\[
 N=x_2\partial_{x_1}+x_3\partial_{x_2}+h\partial_{x_3}.
\]

これは Jordan 型 (4,1) の nilpotent 線形場である。
C では ℓ=w−3n9h、D では ℓ=w−2n9h を、もう一つの超平面の式とする。
両者で n4=n5=n8=0 だから、全係数から

\[
 P^\sharp d\ell=-r\ell N,
 \qquad \Pi^\sharp d\ell=-(r+1/5)\ell N,
 \qquad \Pi=P-\tfrac15E\wedge N
\]

が成立する。従って ℓ を無限遠へ替えると、その affine modular 場は (5r+1)N となる。
ℓ の w 係数は1なので、ker N 上で ℓ=1 となる点を基点に選べる。
これにより定数部は零、線形部は rank三 nilpotent Jordan4 となる。
`MULTIPLE_HYPERPLANE_NILPOTENT_SUPPORT.md` の全支持定理が、全パラメータで適用できる。

さらに C の実因子は `2{h=0}+3{ℓ=0}`。
唯一の三重超平面 ℓ は `im N*` に入らないので、Jordan4 の s=1/5 成分に入る。
D の実因子は `3{h=0}+2{ℓ=0}` で、三重超平面 h は `im(N*)³` に入るから、
Jordan4 の s=1/15 成分に入る。
これは非線形 affine 自己同型を射影座標変換の代わりに用いる議論ではなく、
実際の一次超平面の変更と Euler 鉛直補正による全係数の帰着である。

## 3. EとFの全源は二つの22次元全成分になる

各源の全五自由パラメータを1にした実点を取り、全150二次 lift 係数と
全350射影 Jacobi を保った有限計算を行った。

| 源 | rank D の下界 | vertical | 鉛直・軌道・scaleの下界 | 全実源を加える下界 |
|---|---:|---:|---:|---:|
| E | 103 | 24 | 47 | 47 |
| F | 103 | 24 | 47 | 47 |

記録された非零 minor により下界は標数零でも成り立ち、D が全実方向を零化することは
有理係数で厳密に確認した。従って D の kernel 次元は高々47かつ少なくとも47である。
射影全接次元と実既約族の下界はどちらも22となり、点は滑らかで全成分は唯一。
各全源の PGL 飽和閉包はその22次元成分そのものとなる。
両代表の実軌道も22次元なので、それぞれの成分の稠密開軌道を与える。
別担当の標準 Fraction だけによる独立再構成も受理済みである。
`CLASSIFICATION_NIL31_LONG_SOURCES_REVIEW_20261107.json` は producer・SymPy・NumPy を
import せず、全五源から全 raw Jacobi・Pfaffian・modular・正準 lift と全350Dを再構成した。
全 D103、vertical24、orbit47、実源47の QQ 挟み撃ちを1.44秒で確認した。

### 実因子と相違

E の一般因子は `2H+C3`、F は `3H+Q3`。
E の三次式は x について一次で係数 h²、残部は h=0 で z³/3 となるので primitive。
Gauss の補題により C3 は既約である。
F の二次式 `hy−z²/2+n7h²` は rank三で既約である。
非被約五次因子という閉条件は両全源の閉包に成立し、SNC点を持つ MPS 成分とは異なる。
22以外の次元の成分とも異なる。
F は既存三重超平面成分の rank四二次式と異なり、他の22次元成分とは因子次数・重複度で区別される。

E と同じ `2H+C3` 型を持つ候補は、次の intrinsic な cofactor で区別する。
唯一の二重超平面を ℓ とし `Xℓ=Π♯dℓ/ℓ` とする。
この線形場の Jordan 型と、ℓ の `im Xℓ*` への所属は射影同型不変である。

| 成分の源 | cofactorのJordan型 | 二重超平面の一次式 |
|---|---|---|
| nilJordan4 c=0、s=1/3 | (4,1) | im Xℓ* に入らない |
| 本ノートE、nil(3,1)長い鎖top | (4,1) | im(Xℓ*)³ に入る |
| nil(3,1)短い鎖top、r=−1/3 | (3,2) | Jordan型だけで上二者と異なる |

上二者の代表はいずれも22次元の開軌道なので、射影非同型から成分の相違が従う。
短い鎖topの成分についても同じ開軌道の証明がある。
その源と proof は別担当の `MODULAR_NIL31_SHORT_SUPPORT.md` を用いる。

## 4. Aの全源は実Cartan変形でT01へ入る

A源で `w'=w−n9h` とし、さらに `x'=x−n3w'` とする。
係数を直接比較すると、n3 は消え、c=n2+n3n9、a=n0、k=n7 として

\[
 P=aN^3\wedge N^2+cN^3\wedge N+kN^2\wedge N+W\wedge N,       \tag{3}
\]

\[
 N=y\partial_x+z\partial_y+h\partial_z,\quad
 N^2=z\partial_x+h\partial_y,\quad N^3=h\partial_x,
 \quad W=w\partial_w
\]

となる。N²,N³は線形場の行列の二乗・三乗を意味する。
四場は全て可換で、a≠0 が rank四条件である。

四次元 block `(h,z,y,x)` 上の N を、同じ三角行列の対角に
`0,ε,2ε,3ε` を加えた Nε へ変形する。
(3)で N と全 powers を同時に Nε に置き換え、W はそのままにする。
全 ε で四場は可換だから、得られた二次 bivector は厳密 Poisson である。
ε≠0 では四固有値が異なり、Nε,Nε²,Nε³ は Vandermonde により独立となる。
W と合わせた四次元対角場の空間は Euler 場を含まない。
実際、零固有値の座標で全 powers と W は零、Euler は非零だからである。
従って射影 Cartan 四次元を全て張る。

基底 `(Nε³,Nε²,Nε,W)` における交代係数行列の Pfaffian は −a≠0。
よって ε≠0 の一般元は非退化 log-canonical 構造で、実退化因子は五独立超平面の SNC。
ε=0 は元の全 A源なので、全パラメータを保った実変形と閉性により全 A源はT01に入る。
この全係数恒等式、Cartan の四次元性、Pfaffian 非零は別担当が独立に手で確認した。

## 5. Bの全源を三相異超平面と既約二次因子の層へ変形する

### 5.1 全源の正規形

`W=w−n9 h,X=x−n3 W` は行列式一の affine 変更である。
`a=n0,c=n2+n3 n9,f=n5` と置けば、元の全五変数源は非零定数倍を除いて

\[
\begin{aligned}
 P^{xy}&=a h^2+c hz-fhy+fz^2/2,&P^{xz}&=c h^2,\\
 P^{xw}&=w(fh-y),&P^{yz}&=0,\\
 P^{yw}&=-zw,&P^{zw}&=-hw
\end{aligned}                                                    \tag{8}
\]

となり、h は Casimir である。全 rank四条件は `(a,f)≠(0,0)`。
この正規形は一軌道という主張ではなく、全三変数の実源を保持する。
全係数の普遍恒等式は `CLASSIFICATION_NIL31_LONG_B_CANONICAL_20261107.json` で確認した。

別の Poisson 超平面 w を無限遠に取る。四変数 `(h,x,y,z)` の Euler 場 E と

\[
 K=h\partial_x,\qquad N=y\partial_x+z\partial_y+h\partial_z,
 \qquad R=-hy+z^2/2
\]

を用いれば、(8)へ Euler 鉛直項 `E5∧(fK−N)` を加えた w-Casimir lift は

\[
 Q_0=\Pi_0+E\wedge A_0,\qquad A_0=5fK/4-N,
\]

と書ける。四次元体積形式 `vol=dh∧dx∧dy∧dz` に対して、\(\Pi_0\) は次の構成の ε=0 である。

### 5.2 全パラメータの閉一形式による構成

\[
 H=h+\epsilon z,\quad q_\epsilon=aH^2+fR,\qquad
 N_\epsilon=N+\epsilon(2y\partial_y+z\partial_z)
\]

と置く。次の三次係数一形式は全 `(a,c,f,ε)` について多項式である。

\[
 \omega_\epsilon=q_\epsilon(h\,dz-z\,dh)
       +ch(H\,dR-2R\,dH).                                   \tag{9}
\]

`f ε≠0` の開で説明すると

\[
 \frac{\omega_\epsilon}{hHq_\epsilon}
 =\frac1\epsilon d\log(H/h)
       +\frac cf d\log(q_\epsilon/H^2)                        \tag{10}
\]

であり、右辺は閉じている。また Nε(H)=εH、Nε(qε)=2εqε なので、
右辺の Nε による contraction は1である。
従って `L_Nε ωε=3ε ωε`、`tr Nε=3ε` となる。

\[
 \iota_{\Pi_\epsilon}\mathrm{vol}=\tfrac14 d\omega_\epsilon,
 \quad A_\epsilon=-N_\epsilon+\tfrac{3\epsilon}{4}E+\tfrac{5f}{4}K,
 \quad Q_\epsilon=\Pi_\epsilon+E\wedge A_\epsilon             \tag{11}
\]

と定義する。ωε は Euler に contraction 零で重み四だから、
`ι_E(dωε/4)=ωε`。閉じた(10)から `ωε∧dωε=0` となり、さらに `(dωε)²=0` が従う。
従って Πε は divergence 零・rank高々二の strict Poisson tensor である。
上の Lie 微分恒等式は `[Nε,Πε]=0` を与える。
ωε は x に依存せず dx 成分を持たないので、K も Πε を保つ。
E も斉次二次 Πε と可換であり、従って `[Aε,Πε]=0`、`tr Aε=0`。
これにより Qε は strict Poisson となる。

(9),(11)の係数は多項式なので、導出に使った fε≠0 を除いた全境界でも
Jacobi 恒等式は成り立つ。実際、全パラメータの係数として全四 strict Jacobi、
Πε の Jacobi と divergence、`[Aε,Πε]=0` を直接検算した。
ε=0 では (8) と上記 Euler 差を通じて全係数が一致し、選んだ一点だけの変形ではない。

Pfaffian は厳密に

\[
 \operatorname{Pf}(Q_\epsilon)=-hHq_\epsilon.               \tag{12}
\]

P4 に w-Casimir lift すると、その実五次因子は `w h H qε` である。
二次式 qε の `(h,y,z)` における対称行列式は

\[
 -f^2(2a\epsilon^2+f)/8.
\]

従って各 f≠0 の固定源点について、ε の稠密開では qε は rank三の既約二次式である。
ε≠0 では h,H,w は相異なる三超平面なので、その全点は既済の
**三相異超平面＋既約二次因子の全支持定理**に入る。
その定理の有限個の閉成分の和は閉集合だから、ε=0 の元点も入る。
最後に f≠0 は全 B 源の稠密開であり、全源自体が既約なので、f=0 の rank四境界も同じ有限閉和に入る。
`a=f=0` の rank二境界は既存の全 rank二分類に属する。
これで B の全 rank四支持を被覆し、long-top 全六源の帰属が閉じた。

この証明には、変形全体が一つの既知成分に入るという推定も、元点の局所唯一性も必要ない。
全3H+Q層の有限閉被覆と、各固定源点の実変形を使っている。
証明書は `CLASSIFICATION_NIL31_LONG_B_DEFORMATION_20261107.json`、
producer は `classification_nil31_long_b_deformation.py`、有限契約は
`jobs/nil31_long_b_deformation.json`。実測0.350秒。例外層担当が閉一形式の恒等式、Lie微分、rank二canonical lift、
strict Poisson 構成、全3H+Q定理の適用、ε=0とf=0への閉性を独立に手確認して受理した。
generic-structure担当も閉一形式構成と全支持適用を別に手確認した。
元の全源とのEuler差は公開された全パラメータ恒等式を入力とする。

### 5.3 局所診断の限界と今回の幾何への橋

旧指定点は D rank100、軌道・鉛直・scale の接像46、実源接像47であり、
全二次余核障害は三法線上で零だった。
`CLASSIFICATION_NIL31_LONG_B_NORMAL_QUADRICS_20261107.json` はこの限定的事実を保持する。
また全閉 Q51/Q52 グラフはともに不所属で、
`CLASSIFICATION_NIL31_LONG_B_GRAPH_20261107.json` に全残差がある。
これらを元点の成分性や全支持排除へ拡大していない。

w-scaling 不変な全正規形 (8) では、実源方向を含む真の切片の三法線は全て重み零であった。
`CLASSIFICATION_NIL31_LONG_B_WEIGHTED_SLICE_20261107.json` はこの結果だけを保存し、
重みによる有限な陰関数次数上界は主張しない。
この固定部分を w∞ の斉次二次 tensor として読み直し、(9)の閉有理一形式へ進んだことが
今回の全支持被覆に必要な幾何的変更である。

## 6. 全long-top型のβ非優越はBの帰属に依存しない

固定 affine 場 `A x+c`、A の型 (3,1)、c が長い鎖の top である場合、
Aff4 における安定化群の次元は6である。
実際 GL4 の A 中心化群は6次元で、c mod im A の固定が二条件を課す。
その後、平行移動 v は `Av=Mc` の解として ker A の二次元を持つ。
従って場の軌道は `20−6=14` 次元。
全固定場支持は上の六つの五次元源で尽くされるから、全 affine 支持は高々19次元である。

基点の値が symplectic である開上、GL4 は六次元の非退化定数 bivector 空間に推移的。
従って固定した値 J の fiber は高々 `19−6=13` 次元で、primitive β の像も高々13。
重複超平面を無限遠へ送る perspective が primitive β を保つ既存の橋を用いれば、
この全 modular 型は W=A16 を優越しない。
この非優越の証明は§5の新たな全支持被覆とは独立であり、形式切片の零次元性を仮定しない。

## 有限検算の入口

`code/data/classification_nil31_long_sources.json` は未被覆だった A,B,E,F の全源と各一代表を保持する。
共通 runner `code/classification_modular_source_components.py` は全 Jacobi、全 source、Pfaffian、
全周囲微分を扱い、`CLASSIFICATION_NIL31_LONG_SOURCES_20261107.json` に保存した。
E_multiple_hyperplane の90秒・1536MiB契約 `jobs/nil31_long_sources.json` で実行、四源合計約1.66秒。
すでに幾何的に被覆した C,D に対する接空間の再計算は行っていない。

```sh
python3 scripts/classification.py gate research/classification/jobs/nil31_long_sources.json
ulimit -v 1572864
timeout 90s .venv/bin/python code/classification_modular_source_components.py \
  code/data/classification_nil31_long_sources.json /tmp/nil31-long-sources.json
```
