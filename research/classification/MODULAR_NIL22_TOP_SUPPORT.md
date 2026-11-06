# Jordan (2,2) と長い鎖の定数場：全八源と Cartan 被覆

2026-11-07。E_multiple_hyperplane。全 affine modular 場は
`Z=y∂x+w∂z+∂y`。その斉次化は Jordan 型(3,2)で、h は長さ三の鎖の top covector。
本ノートは全八源の全支持を既知閉包または二つの新成分へ被覆する。一源は T01 に入る。A=−1 は§7の有限な全横断支持判定により
新しい22次元全成分を与える。A=−1/3 も別の22次元全成分となり、
両件は独立 Fraction 実装による全周囲再構成まで受理済み。
A=1/4、A=−2、A=−1/2の主枝と追加枝は実 Poisson 曲線から既知 MPS へ接続する。
A=1 も§11の全パラメータ Aff(3,2) 極限により T19 に入る。従ってこの affine 型の未被覆源はない。

## 1. 全15変数の残式と有限の尽くし

全90係数の線形 reduction は
`CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json` の第五 chart。
自由変数は n0,…,n14。`A=n11,B=n14,u=n10,v=n13` と置く。
全 Jacobi は次の八式と正確に同値である。

\[
\begin{split}
 &(A-B)(2A+1),\
 &(5A-6B+1)u+(A-1)v,\
 &3(v-u)n_9+(4A+2)n_5+3(A-B)n_8,\
 &3un_5+(A+2)n_4+3(A-B)n_7+3n_{12}n_9,\
 &(3A+1)n_3-(3B+1)n_6-3n_5n_9,\
 &(2B+1)n_0-2un_3-(2A+1)n_2+2vn_6+2n_4n_9+2n_5n_8,\
 &vn_0-(A+1)n_1-un_2+n_{12}n_6+n_4n_8+n_5n_7,\
 &3u(u-v)+(4A-1)n_{12}.                                      \tag{1}
\end{split}
\]

A≠−1/2 なら B=A。さらに A≠1 では v=u,n5=0 となり、

\[
 (4A-1)n_{12}=0,\quad (A+2)n_4+3n_{12}n_9=0,
 \quad (3A+1)(n_3-n_6)=0.
\]

これに残りを代入すると、一般 rank4 のためには
`A∈{−1/2,1,1/4,−2,−1/3,−1}` が必要。
この集合の外では n12=n4=n5=n1=0,n3=n6,n0=n2,v=u となり、全 Pfaffian が恒等的に零。

A=−1/2 では (1) を順に解き、

\[
\begin{gathered}
 v=-(1+4B)u,\quad n_{12}=2(1+2B)u^2,\
 n_3=-2(3B+1)n_6-6n_5n_9,\
 n_4=-2un_5+(2B+1)(n_7-4u^2n_9),\
 (2B+1)(n_8+4un_9)=0,\
 (2B+1)(n_0+2un_6+2n_7n_9-8u^2n_9^2)+2n_5(n_8+4un_9)=0.
\end{gathered}
\]

n1 は残る一次式から一意に定まる。
従って B≠−1/2 の主枝と、B=−1/2 で n5=0 または n8=−4un9 の二枝、計三源となる。
正確な全代入は `code/data/classification_nil22_top_sources.json` に保存した。
全八源は七自由変数の既約 affine 源の rank4 開であり、全係数・全残式・raw Jacobi を保持する。

## 2. 全源の射程と接空間検査

| 源 | 全接上界 | 実族下界 | 状態 |
|---|---:|---:|---|
| A=−1/2 主枝 B≠−1/2 | 23 | 22 | §9の普遍 pencil・SNC極限で全 MPS |
| A=B=−1/2,n5=0 | 27 | 22 | §3で全 T01 被覆 |
| A=B=−1/2,n8=−4un9 | 25 | 22 | 別ノートの全分裂 pencil で MPS |
| A=1 | 24 | 22 | 全パラメータ Aff(3,2) 極限で T19 |
| A=1/4 | 23 | 22 | 別ノートの実曲線と SNC 極限で MPS |
| A=−2 | 24 | 22 | §8の実曲線と smoothQ+C3 被覆で MPS |
| A=−1/3 | 24 | 22 | 別ノートで独立検証済み22次元全成分 |
| A=−1 | 23 | 22 | §7で独立検証済み22次元全成分 |

各源の指定 witness は自由変数を JSON の順に1,…,7としたもの一つだけ。
接空間は全350×150 D、実族は全 PGL/vertical/scale と全七 source 方向。
全源の Jacobi/Pfaffian、`DF=0` を正確に確認し、非零 minor は標数零階数下界にだけ用いた。
この表は全局所次元を等号で断定せず、これだけを根拠に成分追加しない。

## 3. A=B=−1/2,n5=0 の全点は T01 に入る

この源の Pfaffian は

\[
 -\tfrac12(n_0-n_2)(w-2u)^2,
\]

なので rank4 は n0≠n2 と同値。
`w'=w−2u,z'=z−2uy` は Z を保つ affine 変更であり、u=0 の表示へ戻す。
その後の全六式は、六つの定数 a,b,c,d,e,f を使って

\[
\begin{split}
 P^{xy}&=a+bw+c(z-yw),\
 P^{xz}&=-yz/2+dw+bw^2+ey+fyw+czw,\
 P^{xw}&=-yw/2,\
 P^{yz}&=-z/2+e+fw+cw^2,\
 P^{yw}&=-w/2,\
 P^{zw}&=-w^2/2.                                             \tag{2}
\end{split}
\]

ここで a=n0+2un6、d=n2+2un6 なので a−d=n0−n2≠0。
Z の中心化群に属する affine 変更

\[
 X=x-2cz+(-2b-4cf)w,\quad Y=y-2cw,\quad
 T=z-2e,\quad W=w
\]

により (2) は

\[
 P^{XY}=a',\quad P^{XT}=-YT/2+d'W+fYW,
 \quad P^{XW}=-YW/2,\quad P^{YT}=-T/2+fW,
 \quad P^{YW}=-W/2,\quad P^{TW}=-W^2/2,
\]

となる。ここで a'=a+2ce,d'=d+2ce、特に a'−d'≠0。

斉次座標 H を戻して、可換な四つの線形場を

\[
 N_3=Y\partial_X+H\partial_Y,\quad N_3^2=H\partial_X,
 \quad N_2=W\partial_T,\quad E_2=T\partial_T+W\partial_W
\]

と置く。冪記号は係数行列の冪から得る線形場を意味する。
全テンソルは正確に

\[
 P=-a'N_3\wedge N_3^2+fN_3\wedge N_2+d'N_3^2\wedge N_2
       -\tfrac12(N_3\wedge E_2+N_2\wedge E_2).                \tag{3}
\]

係数 skew matrix の Pfaffian は `(a'−d')/2`、非零。
三次元 block の N3 を `N3+εY∂Y+2εX∂X` へ、二次元 block の N2 を
`N2+εW∂W` へ変形し、N3² も変形後の行列の平方にする。
E2 はそのまま。全四場は全 ε で可換なので (3) は実 Poisson 族。
ε≠0 では三次元 block の固有値0,ε,2ε、二次元 block の固有値0,εが相異なる。
四場は scalar 行列を法として独立な Cartan 四次元空間を張る。
同時固有座標で (3) は非退化 toric 構造となり、五独立超平面の SNC 因子を持つ。
従って元の全源は閉 toric 成分 T01 に入る。

この全座標変更の六係数、a'とd'、四場分解と Cartan 変形は別担当が独立に手検算した。
点の接空間が27という上界を持つことを、全源の所属証明の代用にはしていない。

## 4. nil31 short-top の一源への帰結

`MODULAR_NIL31_SHORT_SUPPORT.md` の r=1/3 追加枝は実因子 `2H+3L` を持つ。
L=z を無限遠へ取り直すと canonical cofactor からモジュラー場は
非零 scalar 倍の `y∂x+∂y+h∂w`、すなわち現在の nil22 top 型となる。
新しい実因子は `3H+2L`。

全八源の Pfaffian を比べると、この一般 rank4 因子型は §3 の源にしか現れない。
主枝 B≠−1/2 の Pfaffian は xw² の係数が非零で三次。
A=B=−1/2 の追加枝も n5≠0 なら yw² の係数が非零で三次であり、n5=0 は §3 との交わり。
A=1 は三次、A=1/4 と−2はそれぞれ rank3,rank4 の既約二次因子、
A=−1/3 は三次の一一次式の立方、A=−1 は一次 Pfaffian。
従って h が三重で残りが一一次式の平方という型を持つ点は §3 に属する。
これで nil31 short-top の r=1/3 追加枝も全点 T01 へ被覆された。

同じ方法で short-top r=1/2 は現在の A=−1/3 源へ戻る。
`NIL22_RESONANT_AFFINE32_ORBIT.md` はその全七源を一軌道に移し、
全101変数IFTを含む真の二法線切片の根基が原点と証明した。
独立 Fraction 証明書 `CLASSIFICATION_NIL22_THIRD_WEIGHTED_SLICE_REVIEW_20261107.json` は
全D/G/101可逆minorと全IFT・必要な三残式を再構成した。
二つの二次式の係数行列式は−175/2112、残る四次式は第一法線を零として
`−25v⁴/1728` となる。従って根基は二法線の極大idealであり、
全源は22次元全成分となる。これで short-top r=1/2 も全被覆される。

## 5. A=−1 の全源の線形場正規形

この源は全七変数を保ったまま、一つの rank4 軌道の問題へ帰着する。
まず u=n10 とし `W=w−u,Z=z−uy` と置く。この変更はモジュラー場を保つ。
変更後の係数を

\[
 t=n_1,\quad d=n_2+n_6u,\quad b=n_6,\quad c=n_9,\quad
 e=n_7+n_8u+n_9u^2,\quad f=n_8+2n_9u
\]

と書くと、全六式は

\[
\begin{split}
 P^{xy}&=(y^2-2x)/3+d+bw+c(z-yw),\\
 P^{xz}&=t-2(xw+yz)/3+dw+bw^2+ey+fyw+czw,\\
 P^{xw}&=-yw,\quad
 P^{yz}=-(2z+yw)/3+e+fw+cw^2,\\
 P^{yw}&=-w,\quad P^{zw}=-w^2.
\end{split}                                                     \tag{4}
\]

ここで W,Z を再び w,z と書いた。
\(\delta=-3c/4,j=3f,g=-3e/2\)、
\(e_0=-3(b-\delta f)/5,c_0=-3(d-\delta e)/2\) とし、

\[
 X=x+\delta z+e_0w+c_0,\quad Y=y+\delta w,
 \quad T=z+jw+g,\quad W=w
\]

と置く。行列式1のこの変更もモジュラー場を保ち、(4) は正確に

\[
 \begin{array}{lll}
 P^{XY}=(Y^2-2X)/3,&
 P^{XT}=t-2(XW+YT)/3,&P^{XW}=-YW,\\
 P^{YT}=-(2T+YW)/3,&P^{YW}=-W,&P^{TW}=-W^2
 \end{array}                                                     \tag{5}
\]

へ移る。全七係数に関するこの恒等式を、元の公開全源と直接比較して確認した。
斉次座標 h を戻し、再び `(x,y,z,w)` と書けば

\[
 P=t h^2\partial_x\wedge\partial_z+N\wedge S,\qquad
 N=y\partial_x+h\partial_y+w\partial_z,\qquad
 S=(2x\partial_x+y\partial_y-2z\partial_z-3w\partial_w)/3.       \tag{6}
\]

`U=h∂x,V=h∂z` と置くと、N,U,V は可換で
`[S,N]=−N/3,[S,U]=−2U/3,[S,V]=2V/3`。
従って (6) の Jacobi はこの四次元可解 Lie 代数の恒等式としても確認できる。
Pfaffian は tw。rank4 は t≠0 と同値である。
三次元鎖 `(h,y,x)` と二次元鎖 `(w,z)` の相対的な scalar 変更によって任意の非零 t は1へ移る。
したがって全七変数源の rank4 部は一つの PGL 軌道である。
N と S は非可換なので、§3 の可換 Cartan 変形をこの表示へそのまま流用しない。

## 6. A=−1 の未解決範囲を保った有限証拠

`diag(h,x,y,z,w)` の重み `(0,2,1,−1,−2)` の作用で (5) の斉次化は
bivector weight1 の固有元となる。同じ凍結 witness と既存の真の法線を上の affine
変更で移し、重みごとに軌道方向を除くと、3P に対する唯一の法線 weight0 を

\[
 R^{hz}=-yw,\qquad R^{yz}=4h^2-xw,\qquad R^{zw}=4zw
\]

と取れる。全 D の階数102、軌道・vertical・scale の階数47、R を加えた階数48を QQ で確認した。
`D R=0` も全350式で成立する。しかし R の自己 Jacobi は零でなく、
weight−1 の一つの二次補正を用いた `3P+τR+τ²Q` も全 Poisson 族にはならない。
従ってこの計算は新成分性や既知成分所属を証明しない。
同じ法線について追加の高次補正や新しい witness の探索は実行していない。

証拠は `reports/CLASSIFICATION_NIL22_EQUIVARIANT_20261107.json`、
producer は `code/classification_nil22_equivariant.py`、
有限契約は `jobs/nil22_equivariant_deformation.json`。全源の座標恒等式と
同一 witness の equivariant 構成を合わせた実測は1.08秒以内。

候補の除外にも注意を要する。旧 `legacy/v08/sections/cplus_component.tex` の
CPLUS_EXT24 の一般全接次元は25であるから、接次元の上半連続性により、その全閉包の
全点は全接次元25以上である。CMINUS_MAG23 も一般全接次元25で、一般軌道次元は21。
従って上表の全接上界23または24の witness はいずれの成分にも属さない。
軌道22の witness には CMINUS の軌道次元による除外も使える。
因子の次数分割だけでこれらの既知成分への帰属を結論してはならない。

## 7. 有限な全横断支持判定による A=−1 の22次元成分性

§6の不成立は、選んだ二次多項式曲線が Poisson でないという限定的結論であった。
ここでは重みによって**全形式横断切片そのものが有限の多項式問題になる**ことを使う。
任意に jet の精度を増やす仮定は置かない。

### 全形式切片の有限化

150次元二次 lift 空間に、対角作用の後で全体を λ^−1 倍する相対作用を入れる。
これは P を固定し、全 Poisson スキームを保つ。
24次元の Euler 鉛直方向、PGL 作用、scale を合わせた実際の接空間 G は次元47で、
このトーラス作用の安定部分表現である。したがって安定な103次元の線形補空間 S を選べる。
実際の群作用と S の積から全150次元空間への写像は P で微分同型となるので、
S は真の形式横断切片を与える。

S 上の Jacobi 微分の階数は102、核は相対重み−1の一方向 R である。
102個の座標 w と102個の方程式を重みを保って選び、可逆 minor で形式的に消去すると、
`w=φ(τ)` は一意である。一意性から φ は相対作用に関して同変である。
全150座標の相対重みは−8から6なので、各座標成分 φ_i は

\[
 \varphi_i(\tau)=c_i\tau^{-\mathrm{wt}(w_i)}
\]

という次数8以下の一つの単項式、または零である。
重みが非負なら原点条件から零であり、一次係数も `D R=0` により零である。
これにより全形式 IFT 写像の有限な全候補が尽くされ、残る全方程式も有限多項式になる。
この論理は別担当が独立に確認した。

### 完全な一法線切片の非零方程式

`classification_modular_weighted_slice.py` は全150方向から実 G を再構成し、
G の重み安定性、安定補空間103、制限微分102、唯一法線重み−1を QQ で確認する。
実際に用いた全補空間の座標番号、102方程式、可逆 minor は証明書に保存した。
正規化した法線は§6の R の1/4倍である。

一意な102変数消去の後、最初の非零残差に対して、全350行の余核 functional

\[
 \ell=\tfrac{321}{80}e_{15}^*+\tfrac{351}{40}e_{39}^*
       +\tfrac{39}{40}e_{86}^*+e_{228}^*
\]

は全150微分列を零にし、非零値 `455/497664` を取る。
この残差より前の全係数は350式全てで零である。
従って真の一法線形式切片の ideal は非零で、被約支持は原点だけである。
必要な非零残差が得られた時点で計算を止め、正確な冪零厚さは主張しない。

真の実軌道・鉛直・scale の47次元下界と、零次元切片による上界47が一致する。
鉛直24次元と射影scale一方向を除き、P の軌道閉包は22次元の全 Poisson 成分となる。
§5で全七変数源の rank4 部はこの一軌道と証明したので、全源の帰属も閉じる。

実因子は `4H+L`、H と L は相異なる。
従前70の22次元成分の一般因子は、`H+2Q`、`3H+Q`、`2H+C3`、
`2H+L+Q`、`3H+2L` のいずれかであり、今回の成分と異なる。
次元の違う成分、一般 rank二成分とも一致しない。
将来同じ4H+L型を比較するため、唯一の四重 H の canonical cofactor は
Jordan型(3,2)の N の非零倍で、h は `im(N*)²` に入ることも記録する。

永続証明書は `reports/CLASSIFICATION_MODULAR_WEIGHTED_SLICE_20261107.json`、
有限契約は `jobs/modular_weighted_slice.json`。全周囲・全安定切片の実測は1.86秒。
独立検証は `NIL22_WEIGHTED_COMPONENT_REVIEW.md` と
`CLASSIFICATION_NIL22_WEIGHTED_SLICE_REVIEW_20261107.json`。
標準 Fraction の別算術で全D102、実G47、103切片、102可逆minor、全低次係数の零化と
四項functional・非零pairingを再構成した。台帳変更は統合担当が行う。

## 8. A=−2 全源の smoothQ+C3 への実変形

この全七源も既約源の一般点だけでなく、rank4 部の全点が一軌道となる。
u=n10 として、まず `W=w−u/2,Z=z−uy/2` と置く。
変更後の六係数を

\[
 t=n_4,\quad b=n_6,\quad c=n_9,\quad
 d=n_2+bn_{10}/2,\quad
 e=n_7+n_{10}n_8/2+cn_{10}^2/4,\quad f=n_8+cn_{10}
\]

と書く。\(\delta=-c/3,e_0=-(b-\delta f)/4\)、
\(c_0=-(d+tc/3+ce/3)/2\) として、

\[
 X=x+\delta Z+e_0W+c_0,\quad Y=y+\delta W,
 \quad T=Z+fW-e,\quad W=w-n_{10}/2
\]

へ移ると、全源は

\[
 P=N\wedge S+tU\wedge V,\qquad
 N=Y\partial_X+h\partial_Y+W\partial_T,
 \quad S=2X\partial_X+Y\partial_Y-T\partial_T-2W\partial_W,
 \quad U=h\partial_X,\quad V=Y\partial_T+h\partial_W            \tag{7}
\]

となる。全括弧での恒等式、行列式1、Pf(P)=t(YW−hT)は
`CLASSIFICATION_NIL22_MINUS_TWO_CANONICAL_20261107.json` に保存した。
同証明書の整数係数表示は3Pである。
N,U,V は可換、S に関する U,V の重みは−2,+2なので (7) は Poisson。
t≠0 は正確な rank4 条件で、二次元 block `(T,W)` の scalar 変更が全ての非零 t を移す。

古い凍結 witness の像 t=2 において、相対重み `(0,2,1,0,−1)` の
全二法線切片を有限 IFT で求めた。全101選択式は恒等零で、残る ideal は

\[
 (u^2-16v^2/9,\ uv+4v^2/3)=(u+4v/3)(u,v).
\]

と計算される。全IFTグラフに u=−4s,v=3s を代入して得る曲線 P_s は、
全350式を全 s で満たす**実際の**射影 Poisson 曲線である。
被覆には全 residual ideal の独立再計算を要せず、この曲線の小さい恒等式を使えば足りる。

s=1 の実五次因子は非零定数を除いて QC。座標を再び h,x,y,z,w と書くと

\[
 Q=156hx+54hy+162hz+x^2-3xy+48y^2+54xw-162yw,
\]

\[
\begin{split}
 C={}&54w(18hx+x^2-9y^2)+g(h,x,y),\\
 g={}&972h^3+648h^2x+972h^2y+48hx^2+162hy^2
       +x^3+39xy^2+144y^3.
\end{split}
\]

Q は rank5 である。実際、h,z の非退化二次 block を分離した残りの
x,y,w の対称行列は行列式−34992を持つ。
C の w 係数は rank3 の既約二次式であり、h=0,x=3,y=1でこの係数は零、
g は288となる。従って係数と定数項は互いに素で、Gauss の補題により C は
C 上既約な w 一次多項式である。
よって `SMOOTH_QUADRIC_CUBIC_SUPPORT.md` の全支持定理が P_1 に適用される。

全曲線の各係数は上の対角作用の正確な重みを持ち、s≠0 の全点は P_1 と射影同値である。
したがって全て既知 MPS 閉包の和に入り、s→0 の極限 P_0 も入る。
全七源は P_0 の一軌道なので、この源全体の被覆が閉じる。
ここで新しい23次元成分を追加したとは主張しない。

全有限グラフの入力は `CLASSIFICATION_NIL22_MINUS_TWO_WEIGHTED_SLICE_20261107.json`。
実曲線・全Jacobi・因子・smooth/既約条件の小さい永続証明書は
`CLASSIFICATION_NIL22_MINUS_TWO_MPS_20261107.json`。
別担当が全350曲線恒等式、全重み、端点、Qの行列式229582512、
Cのprimitive一次式の条件を直接確認した（実測0.491秒）。
二法線Gröbner計算の正しさを、この被覆の追加前提には用いない。

## 9. A=−1/2 主枝の全 B を MPS へ覆う普遍 pencil

まず u=n10、B=n14 とし、`W=w−2u,Z=z−2uy` に移る。
変更後の係数を

\[
 b=n_6,\ c=n_9,\ e=n_7-4u^2c,\ s=n_5,\qquad
 d=n_2-(24B+10)ub-24usc
\]

と書く。B≠0 上で

\[
 \delta=c/B,\quad j=-2s/(2B+1),\quad e_0=-2(b+\delta s),\qquad
 c_0=\{d+js-2(B+1)e\delta\}/(2B+1)
\]

と置き、`X=x+δZ+e0W+c0,Y=y+δW,T=Z+jW−2e` へ移る。
これは行列式一の affine 変更であり、全七源の整数係数二ベクトルは

\[
\begin{split}
 P_B^{xy}&=0,\\
 P_B^{xz}&=2(2B+1)xw-2(B+1)yz,\\
 P_B^{xw}&=-(2B+1)hz+2Byw,\\
 P_B^{yz}&=-hz,\quad P_B^{yw}=-hw,\quad P_B^{zw}=2Bw^2
\end{split}                                                       \tag{8}
\]

に移る。全係数の正確な恒等式は
`CLASSIFICATION_NIL22_HALF_MAIN_CANONICAL_20261107.json`。
ここでは表示を簡潔にするため、新座標を再び `(h,x,y,z,w)` と書いた。

次の R_B を使う。

\[
 R_B^{hz}=2yz/B,\quad R_B^{hw}=2yw/B,\quad
 R_B^{xw}=-(2B+1)xz/B,\quad R_B^{zw}=z^2.
\]

全350式に対して `P_B+τR_B` は B,τ の全範囲で projective Poisson である。
canonical diagonal weights `(0,2,1,1,0)` では P_B の重みは1、R_B の重みは−1なので、
固定 B で τ≠0 の全点は射影同値となる。

τ=−B/2 を取り、`(h,x,y,z,w)` に重み `(1,0,1,0,0)` を付ける。
全係数の「単項式の重み−二つの微分添字の重み」は非負で、零次極限が存在する。
その極限に行列式非零の座標

\[
 X=x,\quad U=h-y,\quad V=h+y,\quad Z=2w-z,\quad W=2w+z
\]

を入れると logcanonical になる。非零 log 係数は

\[
 \Gamma_{XZ}=-(2B+1),\quad\Gamma_{XW}=2B+1,\quad
 \Gamma_{UZ}=\Gamma_{UW}=1,\quad
 \Gamma_{VZ}=\Gamma_{VW}=-1,\quad\Gamma_{ZW}=-2B.
\]

実 Pfaffian section は `(8B+4)XUVZW` であり、B≠−1/2 なので非退化 SNC。
従って非零 pencil 点は MPS の既知閉包の和に入り、τ→0 で (8) も入る。

B=0 では座標変更に分母が出るが、元の七変数源は B を含めて明示された多項式族である。
その既約源の B≠0 は稠密なので、B=0 の**元源の全像**も上で得た MPS 閉和に入る。
未知 fiber 全体を源閉包から推測したのではなく、同じ既約の実源全体に閉性を適用した。
これで半主枝の全点が被覆される。

全普遍 Jacobi、全 monomial weight、座標行列、toric係数と SNC積は
`CLASSIFICATION_NIL22_HALF_MAIN_MPS_20261107.json` に保存した（実測0.156秒）。
実源 B 方向を含めた有限切片の先行報告は、この被覆の前提ではない。
一般に実源を含む正次元切片から得る局所次元は上界のみであり、共通runnerはこの区別を修正した。
今回に限っては普遍二変数族 `(B,τ)` が真の追加方向を同時に実現するので、
全接上界23と実族下界23が一致することも別途従う。

## 10. 半追加枝の全被覆

`NIL22_HALF_EXTRA_COVERAGE.md` は全七源を a,s の二変数表示に移し、s≠0で
全パラメータの厳密な二次 Poisson pencil を構成した。
その実五次因子は

\[
 -2hw(h+\varepsilon y)
 \{ahw+(a\varepsilon+2s)hz-2syw\}.
\]

`sε(aε+2s)≠0` なら三つの超平面が独立で、残る二次式の階数は4。
したがって既済の全三超平面・既約二次曲面支持定理を適用できる。
ε→0で全元源を戻し、s=0は§3の既済T01源に含まれる。
全 Jacobi/Pfaffian/二次階数の証明書は
`CLASSIFICATION_NIL22_HALF_EXTRA_DEFORMATION_20261107.json`。
別担当による r=z/w 表示の可換場・Lie括弧・閉包適用の独立手検算も受理済みである。
三法線の正次元/零次元判定をこの被覆の追加前提にせず、二次障害計算は実行しなかった。

## 11. A=1 全源の正規形と Aff(3,2) 極限による被覆

u=n10、t=n13−n10 とし、まず `W=w+u,Z=z+uy` と置く。

\[
 b=n_6,\ c=n_9,\ f=n_8-2uc,\ e=n_7-un_8+u^2c,
\]

\[
 d_0=n_2-ub-2tb/3+utc^2/12+tcn_8/3,
 \quad\delta=c/2,\quad e_0=(b-\delta f)/3,
 \quad c_0=(d_0-\delta e)/2,\quad j=-f+\delta t.
\]

行列式一の affine 変更
`X=x+δZ+e0W+c0,Y=y+δW,T=Z+jW` によって、全源は

\[
\begin{array}{lll}
 P^{XY}=2X-Y^2,&P^{XT}=t(X-Y^2/2)+2XW+eY,&P^{XW}=YW,\\
 P^{YT}=e+YW,&P^{YW}=W,&P^{TW}=W(W+t)
\end{array}                                                     \tag{9}
\]

へ移る。全七係数・元の全括弧の恒等式と Pfaffian
`tW(X−Y²/2)` は `CLASSIFICATION_NIL22_ONE_CANONICAL_20261107.json` に保存した。
同証明書の整数係数表示は24Pである。
古い凍結 witness は正確に t=1,e=159 へ移り、新しい点を選んだのではない。

`(X,Y,T,W)→(μ²X,μY,μνT,νW)` と二ベクトルの μ^−1 倍により、
`(t,e)→(νt,μνe)` となる。したがって t e≠0 は一軌道である。
t≠0 が rank4 条件で、e=0は同じ明示的な既約実源の稠密開 e≠0 の閉包に残る。
この正規化の後、以下の実極限が e=0 も含めて全源の所属を証明する。

真の二法線上で、全微分の余核に投影した二次初期空間を一回確認したところ、
独立二次式の ideal は

\[
 (u^2-v^2/101124,\ uv-v^2/318)=(u-v/318)(u,v)
\]

となった。この零次元でない二次初期 ideal だけから全形式支持を推測しない。
`CLASSIFICATION_NIL22_ONE_QUADRATIC_SUPPORT_20261107.json` は全D/実源の階数、
真の法線、全150列を零にする余核 functional と二次式を保存した。
候補直線の任意二次補正は実 Poisson 曲線にならず、三次障害も全二次補正の
曖昧さで消えるため、この局所経路はここで停止した。次の幾何的極限が全源を覆うので、
精密な全形式支持は分類の追加前提としない。

旧 v14 の実族を `Π(b,τ)=N∧H_b+τη_b+τ²θ_b` とする。
旧座標の順序は `(h,y,x,w,z)`、`N=h∂y+y∂x+w∂z`、
`H_b=y∂y+2x∂x+bw∂w+(b+1)z∂z` である。
この**全実族**が23次元 T19 成分へ入ることは、旧成分性証明と SNC 退化の既済入力を使う。
その事実を二つの凍結点の一致から推測したのではない。

\[
 (h,y,x,w,z)=(H,\varepsilon Y,\varepsilon^2X,W,H+\varepsilon Z),
 \qquad b=-1+\varepsilon B                         \tag{10}
\]

とし、座標を移した二ベクトルを ε 倍する。ε≠0で実 PGL・scalar 変換である。
全150係数は ε の多項式となり、ε=0 の極限は

\[
 N_0\wedge H_0+
 \tau\{(Y^2-2HX)\partial_X\wedge\partial_Z
           -3HW\partial_Z\wedge\partial_W\},
\]

\[
 N_0=Y\partial_X+H\partial_Y+W\partial_Z,\qquad
 H_0=2X\partial_X+Y\partial_Y+BH\partial_Z-W\partial_W.
\]

実際、εN∧H_b の極限が第一項、εη_b の極限が波括弧であり、εθ_b は零へ行く。
最後に `Z'=Z+τY` と二ベクトルの符号変更を施すと、正確に (9) の
`t=4τ,e=−B` を得る。任意の t≠0 と e について τ=t/4,B=−e と選べるため、
e=0 を含む**全正規源**が T19 の閉包に入る。

詳しい独立レビュー入口は `NIL22_ONE_AFF32_COVERAGE.md`。
全旧族の projective Jacobi、全150係数の ε 多項式性、極限、最終 shear、
全源恒等式、旧 witness の24倍による移送は
`CLASSIFICATION_NIL22_ONE_AFF32_DEGENERATION_20261107.json` に保存した（0.475秒）。
別担当も N∧H,η,θ の極限と最終六括弧を独立に手展開して受理した。
よって§1の全八源による尽くしと合わせ、この affine modular 型は全点で被覆される。

証拠：`reports/CLASSIFICATION_NIL22_TOP_COMPONENTS_20261107.json`、
共通 runner `code/classification_modular_source_components.py`、
有限契約 `jobs/nil22_top_components.json`。全八源の実測合計5.30秒。
