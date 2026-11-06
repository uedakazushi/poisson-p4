# Jordan (3,1)・短い鎖の定数モジュラー場：全六源と一つの新成分

2026-11-07。E_multiple_hyperplane。座標 `(x,y,z,w)` で
`Z=y∂x+z∂y+∂w` を固定する。全 rank4 支持の全源還元を示し、
三源はそれぞれ22次元全成分に入り、r=1/3追加枝は T01 に入る。
新たな r=1 と r=1/2 の有限支持判定も、独立算術レビューまで受理済みである。
r=1/3主枝と r=2/3 の二源の帰属は未完了。

## 1. 全六源の尽くし

全90係数に対する div/invariance/最高部 wedge 条件の解は14変数 n0,…,n13。
全 Jacobi 残式と全 Pfaffian は
`CLASSIFICATION_MULTIPLE_HYPERPLANE_RANK2_NIL_SOURCES_20261107.json` の第三 chart に保存した。
`r=n13,v=n9` とおく。

残式はまず `n8=0,(3r−1)(r+v)=0,(3r+1)n7=0` を与える。
r≠1/3,−1/3 なら v=−r,n7=0、次いで n2=0 と

\[
 (3r-2)n_3=0,\quad (3r-1)(n_1-n_6)=0,\quad
 (1-2r)(n_4-n_5)=0,\quad (r-1)n_0=0.
\]

r が `{−1/3,1/3,1/2,2/3,1}` に入らなければ、これらを全 Pfaffian に代入して零となる。
したがって rank4 の全点はこの五値に限られる。
r=1/3 の全残式はさらに
`(1+3v)(n6−3n12n10)=0` を与え、二つの既約源に分かれる。

全六源の正確な三角代入と全自由変数は
`code/data/classification_nil31_short_sources.json` に記録した。
どの源も六自由変数の affine 空間の principal open である。
全残式だけでなく、元の四つの affine Jacobi 多項式を全パラメータで再確認した。

| 源 | 自由変数 | rank4 の正確な条件 |
|---|---|---|
| r=1/3 主枝 | n1,n5,v,n10,n11,n12 | v≠−1/3 |
| r=1/3 追加枝、v=−1/3 | n1,n5,n6,n10,n11,n12 | n1≠n6 |
| r=−1/3 | n5,n6,n7,n10,n11,n12 | n7≠0 |
| r=2/3 | n3,n5,n6,n10,n11,n12 | n3≠0 |
| r=1/2 | n4,n5,n6,n10,n11,n12 | n4≠n5 |
| r=1 | n0,n5,n6,n10,n11,n12 | n0≠0 |

個々の代入で原点や一般点だけを選んだのではなく、全源をそのまま保持した。
特に r=1/3 の追加枝を主枝の閉包で自動的に覆ったとは主張しない。

## 2. 全周囲の有限検査

各源で、上表の順に自由変数を1,2,3,4,5,6とする一つの固定 witness を使った。
全元方程式、全150方向、全 PGL/vertical/scale と全六 source 方向を保持する。

| 源 | rank D | rank V | rank M | rank F | 全接上界 | 実族下界 |
|---|---:|---:|---:|---:|---:|---:|
| r=1/3 主枝 | 101 | 24 | 46 | 47 | 24 | 22 |
| r=1/3 追加枝 | 98 | 24 | 45 | 47 | 27 | 22 |
| r=−1/3 | 103 | 24 | 47 | 47 | 22 | 22 |
| r=2/3 | 101 | 24 | 47 | 47 | 24 | 22 |
| r=1/2 | 101 | 24 | 47 | 47 | 24 | 22 |
| r=1 | 102 | 24 | 47 | 47 | 23 | 22 |

ここで表示する D,V,M,F の数は、保存された非零 minor による標数零階数下界。
`DF=0` は整数上で直接確認した。r=−1/3 では上下の一致から階数の等号が従う。
それ以外の行について、接上界だけから全局所次元や既知成分所属を結論しない。
同じ源の witness を取り替える点探索は実行しない。

## 3. r=−1/3 は新しい22次元全成分

この源の全代入は

\[
\begin{gathered}
 n_{13}=-1/3,\quad n_9=1/3,\quad n_8=0,\quad
 n_2=-3n_{11}n_7/2,\quad n_3=-n_{12}n_7,\\
 n_1=n_6+3n_{12}^2n_7/2,\quad
 n_4=n_5+3n_{11}n_{12}n_7,\\
 n_0=\tfrac34n_7(n_{10}n_{12}+\tfrac32n_{11}^2+n_6).
\end{gathered}                                                   \tag{1}
\]

自由変数は n5,n6,n7,n10,n11,n12、n7≠0。
表の真の全接上界22と実族下界22から witness は滑らかで唯一の全空間成分を通る。
全六源の PGL 飽和像閉包は既約で22次元なので、その全成分そのものである。
さらに軌道だけでも射影次元22を持つから、この成分は稠密開軌道を持つ。

実 Pfaffian は n7 倍の三次式である。affine 線形座標変更

\[
\begin{split}
 W&=w-n_{12}z-3n_{11}/2,\
 Y&=y,\
 T&=z,\
 X&=x+\tfrac34(n_{10}n_{12}+n_6)z+\tfrac32n_{10}n_{11}+n_5
\end{split}
\]

で、正確に

\[
 \operatorname{Pf}P=n_7(X-YW+TW^2/2).
\]

従って実五次因子は `h²C3`、
`C3=h²X-hYW+TW²/2`。
C3 は X に関して一次で、その係数 h² と残りの項が互いに素なので既約である。
特に重複超平面は h=0 一つだけである。

同じ一般因子型を持つ Jordan4 源の新22次元成分との区別には、因子型だけでは足りない。
canonical divergence-free representative は
`Π=P−E∧Z_hom/5`、

\[
 Z_{\rm hom}=y\partial_x+z\partial_y+h\partial_w.
\]

唯一の二重超平面 h の cofactor は `Π♯dh/h=−Z_hom/5`。
この線形場の Jordan 型は **(3,2)**、冪零指数は3。
Jordan4・定数部零の源では型 **(4,1)**、冪零指数は4である。
もう一つの nil31 long-top 源も型(4,1)である。
唯一の二重超平面から自然に定まる cofactor の冪零指数は PGL 同値と全体の scalar 倍で不変なので、
それぞれの稠密開軌道は一致しない。従って成分も相異なる。

他の同次元非MPS成分とは一般因子型で区別され、次元が異なる成分とは同一でない。
MPS は最小でも23次元。一般 rank2 成分とは rank 条件で区別する。
独立再構成は `reports/CLASSIFICATION_NIL31_SHORT_COMPONENT_REVIEW_20261107.json`。
全パラメータの raw Jacobi、QQ 上の D103/V24/M47/F47、全実方向の D 消滅、
canonical cofactor の冪の階数3,1,0、三次因子の既約性を独立に確認した。
正式な台帳追加は統合担当が行う。

## 4. 他の五源と、残る r=2/3 源

他の五源の一般因子型は順に
`2H+L+Q_rank3`、`2H+3L`、`2H+L+Q_rank4`、`3H+2L`、`4H+L`。
このうち r=1/3追加枝の `2H+3L` は、L を無限遠に取る正確な変更により
`MODULAR_NIL22_TOP_SUPPORT.md` §3 の全 Cartan 源へ戻り、全点が T01 に入る。
同ノート §4 は八源の全 Pfaffian を比較して、この橋の補集合がないことも確認する。
さらに r=1 は§6で全成分へ、r=1/2 は nil22 A=−1/3 の全源への正確な橋と
`NIL22_RESONANT_AFFINE32_ORBIT.md` §4 の有限二法線支持判定により、別の22次元全成分へ入る。
後二件も別担当による独立算術レビューが完了しており、台帳への同期は統合担当が行う。
r=1/3 主枝も§8の全パラメータ変形で被覆する。残る r=2/3 源は§7へ正規化したが、
因子型が既知族と同じという理由だけで所属を主張しない。全 source の実変形または
別の幾何的認識が次の義務である。
接空間の二次以上を際限なく精密化する経路は採らない。

## 5. r=1 全源の一軌道化と有限な安定トーラス切片

r=1 の六自由変数を n0,n5,n6,n10,n11,n12 とし、

\[
 b=3n_{10},\quad f=3n_{12}/5,\quad g=3n_{11}/2,
 \quad c=(3bn_{12}-3n_6)/4,\quad e=-3n_5+3bn_{11}
\]

と置く。モジュラー場 `N=y∂x+z∂y+h∂w` を保つ行列式一の affine 変更

\[
 X=x+cz+bw+e,\quad Y=y+b,\quad Z=z,\quad W=w+fz+g
\]

により、**全六パラメータ源**の全括弧は

\[
\begin{array}{lll}
 P^{XY}=t-(XZ+Y^2)/3,&P^{XZ}=-YZ,&P^{XW}=-X/3+2YW/3,\\
 P^{YZ}=-Z^2,&P^{YW}=Y/3+2ZW/3,&P^{ZW}=Z,
\end{array}\qquad t=n_0.                                     \tag{2}
\]

へ移る。全元係数とのこの恒等式、raw Jacobi、Pf(P)=tZ、行列式1を
`CLASSIFICATION_NIL31_SHORT_CANONICAL_20261107.json` に記録した。
本節の P は、同証明書の整数係数 `canonical_brackets` に表示した二ベクトルの1/3倍である。
従って証明書の括弧の Pfaffian は9倍、射影点と全支持の議論は同じである。
旧 witness の像が正確に t=1 の (2) になることも確認しており、新しい点の探索ではない。

斉次化した (2) は、線形場

\[
 N=Y\partial_X+Z\partial_Y+h\partial_W,\qquad
 S=(X\partial_X-Y\partial_Y-3Z\partial_Z+2W\partial_W)/3
\]

を使って `P=N∧S+t h²∂X∧∂Y` と書ける。
`[S,N]=−2N/3` で、残る二場 `h∂X,h∂Y` の S 重みは互いに反対。
全 rank4 条件は t≠0 であり、三次元 Jordan block `(X,Y,Z)` を
二次元 block `(h,W)` に対して scalar 倍する作用が t を非零平方倍する。
C 上では全ての非零 t が移り合い、全六源は一つの PGL 軌道となる。

この点は対角重み `(h,X,Y,Z,W)=(0,1,0,−1,1)` に関して bivector weight1である。
P を固定する相対作用の全150座標の重みは−4から3となる。
唯一の真の法線の相対重みは−1と QQ で確認され、
`MODULAR_NIL22_TOP_SUPPORT.md` §7 の有限 IFT 切片の論理が適用できる。

## 6. r=1 の22次元全成分と、nil22 成分との相違

全150係数から、実軌道・Euler鉛直・scale の47次元接空間を再構成する。
そのトーラス安定補空間103次元において、全 Jacobi 微分は階数102である。
唯一の法線 τ の相対重み−1と全周囲の重み範囲から、一意な102変数 IFT 写像は
次数四以下の多項式に尽くされる。
この全多項式を求めて代入すると、選んだ102式は**恒等的に**零であり、
残る全350式の ideal の根基は `(τ)` となる。
実際、保存した残式中に非零定数倍の τ³ があり、全残式は原点で零である。
これは任意精度の jet からの推定ではなく、全形式横断切片の有限な厳密表示である。

真の横断被約支持が原点なので、全局所次元は実軌道等の47次元に等しい。
Euler鉛直24次元と射影scale一方向を除くと、P の軌道閉包は22次元の全成分になる。
§5で全六源の rank4 部が一軌道と確認したので、全源もこの成分で覆われる。
全 nilpotent thickness の分類は主張しない。

一般実因子は `4H+L`。同じ因子型を持つ
`MODULAR_NIL22_TOP_SUPPORT.md` §7 の新22次元成分とは、唯一の四重超平面の
canonical cofactor とその双対像で区別する。現在の場は

\[
 N=Y\partial_X+Z\partial_Y+h\partial_W,
 \qquad (N^*)^2(V^*)=\mathbb C Z.
\]

従って H の式 h は `im(N*)²` に入らない。
nil22 型の場 `Y∂X+h∂Y+W∂T` では、逆に h が `im(N*)²` を張る。
どちらも Jordan 型(3,2)なので型だけでは足りないが、この flag 条件は
PGL 変換と cofactor の非零 scalar 倍に不変である。
両者は22次元の稠密開軌道を持つため、軌道の非同値から全成分の相違が従う。
それ以外の従前70成分とは、22次元に同じ一般因子型がないことと次元・一般階数で区別する。

永続証明書は `CLASSIFICATION_NIL31_SHORT_WEIGHTED_SLICE_20261107.json`。
共通 producer `classification_modular_weighted_slice.py` の入力と全 bytes、
103切片・102可逆minor・法線・全IFT多項式・全残idealを保存した。
有限契約は `jobs/nil31_short_weighted_slice.json`、実測1.45秒。
独立検証 `CLASSIFICATION_NIL31_SHORT_WEIGHTED_SLICE_REVIEW_20261107.json` は
全六源の座標正規化、全D/G/IFTと全350残差を別に再構成した。
特に行57の残式 `440τ³/28431` が非零である。実測2.495秒で受理済み。
台帳への同期は統合担当が行う。

検証資料は `code/classification_modular_source_components.py`、上記六源 JSON、
`reports/CLASSIFICATION_NIL31_SHORT_COMPONENTS_20261107.json`。
有限契約は `jobs/nil31_short_components.json`。
全六源の raw Jacobi/Pfaffian と一つずつの全接空間検査の実測合計は3.54秒。
この report は全 nilrank2 の被覆を主張しない。


## 7. r=2/3 全源の正規化と未決の真二法線

この節は全支持の証明ではなく、残る六変数源を正確に二変数へ帰着する。
元の自由変数 n3,n5,n6,n10,n11,n12 に対し

\[
\begin{gathered}
 t=n_3,\qquad d=n_5-3n_{10}n_{11},\\
 X=x-\tfrac32\{n_6+n_{12}(n_3-3n_{10})\}z+3n_{10}w,\\
 Y=y+3n_{10},\quad Z=z,\quad W=w+n_{12}z+3n_{11}
\end{gathered}                                                    \tag{14}
\]

と置く。行列式一の affine 変更であり、全六括弧は

\[
\begin{aligned}
 P^{XY}&=dZ-tY+2tZW-Y^2/3,&P^{XZ}&=tZ-2YZ/3,\\
 P^{XW}&=d+YW/3,&P^{YZ}&=-2Z^2/3,\\
 P^{YW}&=(Y+ZW)/3,&P^{ZW}&=2Z/3.                  \tag{15}
\end{aligned}
\]

に一致する。元証明書は整数倍 3P を記録する。
全パラメータの係数恒等式、全 affine Jacobi、行列式と既存 witness の移送は
`CLASSIFICATION_NIL31_SHORT_TWO_THIRDS_CANONICAL_20261107.json` に保存した。
実五次因子は非零定数倍の `t h² Z(hY−ZW)` で、t≠0 が rank四の条件である。

座標を `(h,X,Y,Z,W)↦(h,μ²νX,μνY,νZ,μW)` とし bivector を μ⁻¹ 倍すると
`(t,d)↦(μνt,μ²νd)` となる。したがって td≠0 は一つの PGL・scalar 軌道である。
d=0 は同じ明示的既約多項式源の閉集合として残す。
一般点の既知閉包への所属が得られた場合だけ、閉性でこの集合へ戻すことができる。

元の witness は t=1,d=−58 へ送られる。この同じ点で全350×150微分 D、
軌道・vertical・scalar の G、および二つの実源方向を保持した真の切片を構成した。
QQ で `rank D=101,rank G=47` で、二法線の全二次余核像は

\[
 (u^2-v^2/84100,\;uv-v^2/290)
       =(u-v/290)(u,v).                              \tag{16}
\]

これは初期支持に一本の直線を残す。全形式支持が直線という結論や成分性は含意しない。
`CLASSIFICATION_NIL31_SHORT_TWO_THIRDS_QUADRATIC_20261107.json` は
全Dの左余核関数、独立二次像、真の法線代表を保存する。
ここから高次 jet を際限なく延長せず、既知族への実変形または幾何的認識を次の義務とする。

幾何を比較するため、開集合 Z≠0 で
`q=Y/Z−W,F=X−YW+ZW²/2` と置くと、括弧は

\[
 \{F,q\}=-2tq,\quad\{F,Z\}=tZ,\quad\{F,W\}=d,
 \quad\{q,W\}=-q/3,\quad\{Z,W\}=2Z/3,\quad\{q,Z\}=0.
\]

これは有理座標による線形化であって PGL による成分所属の証明ではない。


## 8. r=1/3 主枝の全支持は半単純階数二へ変形される

v=n9,u=n10,f=n12,g=n11 とし、
`d=n5−3ug,a=n1−3uf` と置く。v(3v+1)≠0 上で

\[
 X=x+3uw+3d+\frac{a}{3v+1}z,\quad
 Y=y+3u,\quad Z=z,\quad W=w-\frac f v z                 \tag{17}
\]

は行列式一の affine 座標変更である。全六源の括弧は正確に

\[
\begin{aligned}
 P^{XY}&=(2v+1)XZ-(v+1)Y^2/2,&P^{XZ}&=vYZ,\\
 P^{YZ}&=vZ^2,&P^{XW}&=X/3+gY,\\
 P^{YW}&=Y/3+gZ,&P^{ZW}&=Z/3                         \tag{18}
\end{aligned}
\]

となる。整数係数を保存した旧源は6Pの規約である。
`CLASSIFICATION_NIL31_SHORT_THIRD_MAIN_CANONICAL_20261107.json` は、
全六変数を残した座標恒等式・全 Jacobi・行列式・旧 witness の正確な移送を確認する。
以下では X,Y,Z,W を再び x,y,z,w と書く。

三変数 q=2xz−y²、E3=x∂x+y∂y+z∂z に対し

\[
 q_\varepsilon=2xz-y^2+\varepsilon x^2,\qquad
 N_\varepsilon=y\partial_x+(z+\varepsilon x)\partial_y
\]

と置く。三変数二ベクトル π_ε を

\[
 \iota_{\pi_\varepsilon}(dx\wedge dy\wedge dz)
       =\frac12\{vz\,dq_\varepsilon+(v+1)q_\varepsilon\,dz\}
                                                               \tag{19}
\]

で定め、次の実多項式族を使う。

\[
 P_\varepsilon=\pi_\varepsilon+
       (E_3/3+gN_\varepsilon)\wedge\partial_w.                 \tag{20}
\]

(19)の一形式 α は α∧dα=0 なので π_ε は Poisson。
`Nε(qε)=Nε(z)=div Nε=0` より Nε は π_ε を保ち、
π_ε が二次斉次なので E3 も π_ε を保つ。
したがって (20) は任意の v,g,ε に対して真の二次 Poisson 構造で、P0 は (18) である。
その modular 場と実因子は直接

\[
 \operatorname{div}P_\varepsilon=N_\varepsilon+\partial_w,
 \qquad D_\varepsilon\doteq(3v+1)h^2zq_\varepsilon.           \tag{21}
\]

線形部 Aε の特性多項式は `λ²(λ²−ε)` で、ε≠0 なら零固有空間は二次元。
従って Aε は半単純階数二である。定数成分 ∂w は im Aε の外に残るので、
これを零に置き換えず、`MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §2 の
**定数成分を含む全半単純階数二支持定理**を適用する。
3v+1≠0 なら全 ε で rank四であり、generic ε の全点が既知 MPS 閉和へ入る。
ε=0 の閉性により (18) の任意 g も同じ閉和に入る。

(17)で除いた v=0 は、元の六変数**多項式実源**の内部にある。
v≠0 はこの既約源で稠密なので、その全像の閉性により v=0 も被覆される。
3v+1=0 は主枝の局所化から外れ、既に§4で T01 へ被覆した追加枝が保持されている。
この論理は未知 fiber を正規源の閉包から推測するものではない。

別担当も (19)–(21) の全係数と半単純性を手で再導出して独立に受理した。
全パラメータ Jacobi、Pfaffian、modular、特性多項式の小証明書は
`CLASSIFICATION_NIL31_SHORT_THIRD_MAIN_DEFORMATION_20261107.json` に保存した（0.101秒）。
この被覆では局所の接空間・高次障害を追加前提としない。
