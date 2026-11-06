# Jordan (2,2) と長い鎖の定数場：全八源と Cartan 被覆

2026-11-07。E_multiple_hyperplane。全 affine modular 場は
`Z=y∂x+w∂z+∂y`。その斉次化は Jordan 型(3,2)で、h は長さ三の鎖の top covector。
本ノートは全八源を保持し、一つの全源を T01 に被覆する。他の七源の帰属は未完了。

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
| A=−1/2 主枝 B≠−1/2 | 23 | 22 | 未被覆 |
| A=B=−1/2,n5=0 | 27 | 22 | §3で全 T01 被覆 |
| A=B=−1/2,n8=−4un9 | 25 | 22 | 未被覆 |
| A=1 | 24 | 22 | 未被覆 |
| A=1/4 | 23 | 22 | 未被覆 |
| A=−2 | 24 | 22 | 未被覆 |
| A=−1/3 | 24 | 22 | 未被覆 |
| A=−1 | 23 | 22 | 未被覆 |

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

同じ方法で short-top r=1/2 は現在の A=−1/3 源へ戻るが、その帰属はまだ未決。
源への帰着と既知成分への帰属を混同しない。

証拠：`reports/CLASSIFICATION_NIL22_TOP_COMPONENTS_20261107.json`、
共通 runner `code/classification_modular_source_components.py`、
有限契約 `jobs/nil22_top_components.json`。全八源の実測合計5.30秒。
