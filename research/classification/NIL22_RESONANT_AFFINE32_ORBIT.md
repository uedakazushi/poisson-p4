# Jordan (2,2) 定数型 A=−1/3：全源の一軌道化と Aff(3,2) 共鳴

2026-11-07。E_multiple_hyperplane。
`MODULAR_NIL22_TOP_SUPPORT.md` の A=−1/3 全七変数源、同値に
`MODULAR_NIL31_SHORT_SUPPORT.md` の r=1/2 全源を扱う。
以下で全源を一つの明示軌道へ正確に帰着し、全有限横断支持の判定から新しい22次元全成分を得る。
全源の座標変更は別担当が独立に手確認した。全周囲の数値証拠も別の Fraction 実装で独立再構成し受理済みである。

## 1. 全係数を残す中心化群による変更

固定 modular 場を `Z=y∂x+w∂z+∂y` とする。
元の全源の u=n10 に対して `w'=w−3u,z'=z−3uy` は Z を保ち、u=0 の表示へ送る。
この最初の変更は全パラメータで定義され、逆も affine である。
u=0 で残る六係数を a,b,c,d,e,f と書くと、全六括弧は

\[
\begin{aligned}
 P^{xy}&=a+bw+c(z-yw)+(2x-y^2)/9,\\
 P^{xz}&=aw+dy+ew^2+fyw+czw+(2xw-4yz)/9,\\
 P^{xw}&=-yw/3,\\
 P^{yz}&=d+fw+cw^2-4z/9+yw/9,\\
 P^{yw}&=-w/3,\qquad P^{zw}=-w^2/3.
\end{aligned}                                                    \tag{1}
\]

Pf(P)=(e−b)w³/3 であり、元の源で e−b=n3−n6。
したがって一般 rank四の条件は正確に e≠b である。
次の六つの量を置く。

\[
 A=-9c/2,\quad B=9(Af-b),\quad
 s=9(a-Ad)/2,\quad\delta=9d/4.
\]

Z の中心化群に属する、行列式一の affine 変更

\[
 X=x+Az+Bw+s,\quad Y=y+Aw,\quad
 T=z-\delta-9fw,\quad W=w                                 \tag{2}
\]

により、(1)の全六係数は正確に

\[
\begin{aligned}
9P^{XY}&=2X-Y^2,&9P^{XT}&=\lambda W^2+2XW-4YT,\\
9P^{XW}&=-3YW,&9P^{YT}&=-4T+YW,\\
9P^{YW}&=-3W,&9P^{TW}&=-3W^2,
\qquad\lambda=9(e-b)\ne0.                                \tag{3}
\end{aligned}
\]

順に X,Y の shear で c,b を消し、T の平行移動で d、X の平行移動で a、
最後の T の shear で f を消すと (2) を得る。
三回の途中表示でも e−b は保たれる。
全七パラメータの元源から (3) までの六係数恒等式と Pfaffian は
`CLASSIFICATION_NIL22_RESONANCE_NORMALIZATION_20261107.json` で QQ 上確認した。

さらに (W,T) を同時に k 倍すると λ は λ/k となり、ほかの項は変わらない。
したがって λ≠0 の全源は一つの PGL 軌道である。
既存の全周囲証明書で、この軌道は射影次元22、全接空間は次元24である。
この段階の接空間だけから局所次元22や成分性を推論しない。§4で真の横断支持を使って成分性を証明する。

## 2. 共鳴する四次元可解 Lie 代数による表示

斉次座標 h を戻して

\[
 N=Y\partial_X+h\partial_Y+W\partial_T,\quad
 D=2X\partial_X+Y\partial_Y+4T\partial_T+3W\partial_W,
 \quad U=W\partial_X,\quad V=W\partial_T
\]

と置けば、(3)は

\[
 9P=D\wedge N+\lambda U\wedge V.                           \tag{4}
\]

N,U,V は相互に可換で、`[D,N]=−N,[D,U]=U,[D,V]=−V`。
従って (4) は全 λ で厳密に Poisson である。
標準の Aff(3,2) 座標順を (h,Y,X,W,T) とすると、D=H3 は重み (0,1,2,3,4) を持つ。
一般 H_b では `[H_b,U∧V]=(b−3)U∧V` だから、λ≠0 の項は b=3 に限る共鳴拡張である。
全実因子は `2{h=0}+3{W=0}`。他の重複超平面を無限遠とすれば
nil31 short r=1/2 の `3H+2L` 表示へ戻る。

h=1,W≠0 において

\[
 q=X-Y^2/2,\qquad s=T/W-Y,\qquad t=Y
\]

と置くと、(4) は線形 Poisson 括弧

\[
 \{q,t\}=2q,\quad\{s,t\}=s,\quad
 \{W,t\}=3W,\quad\{q,s\}=\lambda W
\]

となる。これは Heisenberg Lie 代数と重み (2,1,3) の導分の半直積である。
この座標は双有理であって PGL ではないので、線形 Poisson の分類から
射影成分への所属を推定してはいけない。

## 3. 共通のrank二境界からは所属を推論しない

v14 の普遍 Aff(3,2) 変形 `N∧H_b+tη_b+t²θ_b` は rank二の境界を
既存 MPS 成分へ入れる。しかし、その証明だけで新項 λU∧V を含む (4) の所属は従わない。
λ=0 の共通境界への退化を、λ≠0 の点の同一成分性と取り違えない。

また、(4) 自体を一般 toric 構造へ変形して T01 に入れる経路は成立しない。
この点の全接空間上界24は、T01 の成分次元25より小さいからである。
線形 Poisson 表示の導分重み (2,1,3) を (a,b,a+b) へ変えるだけでも、
元の括弧 PXT に `(b−a/2)Y²W` という三次項が生じる。
ほかの三次係数を持たないこの項は射影二ベクトルの許容する Euler 型三次尾部ではなく、
そのままでは P4 上の変形にならない。

二つの接方向の存在や一軌道化だけでは、既知成分への包含も新成分性も従わない。
本源では以下の全周囲の支持上界と実族下界が一致する。

入力は `CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json` と
`code/data/classification_nil22_top_sources.json`。
有限検算は `code/classification_nil22_resonance_normalization.py`、契約は
`jobs/nil22_resonance_normalization.json`。全入力 hash と全座標変更は出力に保存した。
これは手導出した普遍恒等式の producer 検算であり、別実装による独立再現とは呼ばない。


## 4. 真の二法線切片の被約支持は原点だけ

(3)を9倍した tensor について λ=1 とする。
これは元の凍結 witness n3=2,n6=3（λ=−9）へ §1 の変更と
(W,T) の同時 −9 倍を適用したもので、新しい標本を選んでいない。
対角 pushforward の座標重みを

\[
 (\operatorname{wt}(h),\operatorname{wt}(X),\operatorname{wt}(Y),
   \operatorname{wt}(T),\operatorname{wt}(W))=(0,2,1,3,2)
\]

とすると P の重みは1。全体を対角パラメータの逆数倍する相対作用は P を固定する。
全150係数で Jacobi 微分 D の階数は101、実 PGL・鉛直・scale の接空間 G の次元は47。
全 G は重み安定であり、103次元の安定線形補空間 S を持つ。
`CLASSIFICATION_NIL22_THIRD_WEIGHTED_SLICE_20261107.json` は G の pivot 行、
全 S 座標、101次の非零 QQ minor と選択方程式を保存する。
真の二法線はいずれも相対重み−1である。

S 上の101変数を一意な形式陰関数で消去し、残る法線を (u,v) とする。
相対 ambient 重みは [−6,4] なので、陰関数の各座標は (u,v) の次数6以下の
斉次多項式または零となる。同変性は形式陰関数の一意性による。
従って必要な全単項式を有限に列挙でき、精度を任意に選んだ有限 jet ではない。
全次数6までの補正を入れた後、選択した101方程式は厳密な多項式として恒等零となった。

残る全350式を代入した ideal の Groebner 基底は

\[
 (v^4,u^2,uv).
\]

被約支持に必要な証拠は、その全 Groebner 計算より小さい。
四次係数の標準並び順を公開 report と同じとすると、行115と157は正確に

\[
 f_{115}=-\frac5{1056}u(35u+36v),\qquad
 f_{157}=\frac5{48}u(7u+12v).                              \tag{5}
\]

この二式の定数係数線形結合から u²,uv が得られる。
行98は四次式

\[
 -\frac5{204578390016}\bigl(
 64917369335u^4+116695455980u^3v+13944847216u^2v^2
 -34710447552uv^3+591951360v^4\bigr).                       \tag{6}
\]

である。(5)を法とした (6) は非零定数倍の v⁴だから、真の切片 ideal の根基は (u,v)。
この結論に exact な非被約厚さの決定は不要である。
独立検証は全 D と真の S、次数4までの一意な IFT 補正、(5),(6) の三式だけでも足りる。
全陰関数の有限化と行の重み同変性により、後の補正がこれらの斉次式へ寄与しないことが分かる。

したがって真の形式切片は被約点である。
47次元の実作用とこの切片の積から全150次元空間への形式同型を使えば、
全局所支持の次元は47。鉛直24次元と射影 scale 一次元を引いて22となる。
§1の22次元軌道閉包は、全 Poisson 空間の22次元既約成分である。
全七変数源の rank四部分がこの一軌道なので、元の全源もこの新成分に入る。
同値な nil31 short r=1/2 全源も被覆された。

実因子は `2H+3L`。新成分の一般点はこの稠密軌道に属する。
同じ因子型の従前22次元 nilJordan4 成分 s=1/5 と s=1/15 は、いずれも唯一の三重 L の canonical cofactor が
Jordan 型(4,1)である。一方、現在の三重 W の cofactor は、規約に合わせて
`−6N/5` であり Jordan 型(3,2)である。
これは PGL と非零 scale を保つ内在的不変量なので、両成分は相異なる。
他の旧22次元成分とは一般因子型が異なり、異なる次元・一般 rank二の成分とも一致しない。

producer は CLI 化した `code/classification_modular_weighted_slice.py`、入力は
`code/data/classification_nil22_third_weighted_slice.json`、有限契約は
`jobs/nil22_third_weighted_slice.json`。全有限判定は6.381秒、入力 hash と全 producer bytes を保存した。
`CLASSIFICATION_NIL22_THIRD_WEIGHTED_SLICE_REVIEW_20261107.json` は別の Fraction 実装で、
全 D101・G47・二法線・真の切片、次数4までの陰関数補正と三残式を独立に再構成し、
根基 (u,v) を受理した。台帳と主原稿への統合は統合担当者が行う。

§1の全六係数、D∧N+λUV 表示、Lie 括弧、rank四条件、block scale による一軌道化は
例外層担当が独立に手で再導出して確認した。初段 u 変更の全源恒等式は公開 producer 証明書への依存である。
