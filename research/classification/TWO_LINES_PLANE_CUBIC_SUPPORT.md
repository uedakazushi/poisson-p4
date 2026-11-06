# 平面三次曲線と二直線：全十三配置・全磁気源と九つの23次元成分

2026-11-07。課題 `E_reduced_nonnormal`。
実退化式が相異なる二超平面と既約三次因子の積で、全五次式の本質変数数が高々三の場合を扱う。
全支持を十三個の明示的な六パラメータ源のPGL飽和へ還元する。
そのうち六源の閉包は真正一・二法線の二次障害から23次元の全Poisson成分と独立検証済みである。
さらに三源について、固定した真の切片の三次初期式から同じ結論を得て、独立算術検証も完了した。
残る四源の成分帰属は未完了であり、型1+1+3全体の被覆完了は主張しない。

## 1. 全係数を保持するfree平面五次曲線への還元

本質変数が二以下の三次式はC上で可約なので、対象の本質変数数は正確に三である。
二超平面をh=0,x=0とし、残る本質変数をy、fiber座標をz,wとする。
h-Casimir代表を取ると実退化式は

\[
 D=h x C(h,x,y),\qquad C\text{ は既約な平面三次式}.             \tag{1}
\]

各因子はPoisson-normalである。h=1上では
Cが{ x,C }=P_xy C_yを割る。CとC_yは互いに素で、P_xyの次数は高々二だからP_xy=0。
従って全構造は

\[
 P=U\wedge\partial_z+V\wedge\partial_w+B\partial_z\wedge\partial_w,
                                                               \tag{2}
\]

となる。初めはU,Vの係数にz,w依存を許す。
Pfaffian条件はdet(U,V)=非零定数×xCである。
z=w=0へ評価すれば、h,x,Cを接線方向に保つ二つの全二次anchorを得る。
Euler場Eと合わせた3×3行列式は非零定数×hxCであり、Saitoの行列式判定から
平面五次曲線hxCはexponents(1,2,2)のfree divisorとなる。

このfree基底をR=C[h,x,y]からR[z,w]へ拡張する。
二次anchorをE,U,Vで書いたときU,Vの係数は全次数零、従って定数である。
h成分が零であることからE係数も零となり、もとのU,Vのfiber依存は消える。
一定のGL2 fiber変更により、選んだ基底U,Vそのものへ戻せる。
この議論はfiber依存を初めから除いた切片だけを扱うものではない。

[U,V]=cU+dVとするとc,dはh,x,yの斉次一次式であり、全Jacobiは
B_z=-c,B_w=-dと同値である。従って全源は正確に

\[
 P=U\wedge\partial_z+V\wedge\partial_w+
    \{-cz-dw+R(h,x,y)\}\partial_z\wedge\partial_w,
 \quad R\in H^0(\mathbb P^2,\mathcal O(2)).                    \tag{3}
\]

Rの全六係数を保持する。逆に(3)は全RでPoissonであり、実因子は(1)のままである。

exponents(2,2)のgradient syzygy二本を用いるHilbert–Burch分解から、
全Tjurina数は

\[
 \tau(hxC)=(5-1)^2-2\cdot2=12.                                \tag{4}
\]

以下の配置分類はこの必要条件を使う。
最後に各残存配置の全二次syzygyと行列式を計算するので、必要条件だけを十分条件と取り違えない。

## 2. 交点数による全配置の有限化

sを三つの大域成分間の相異なる交点数、bを二直線の交点がC上にあるとき1、他を0とする。
局所の積公式μ(fg)=μ(f)+μ(g)+2I(f,g)-1とBézoutから

\[
 \sum_p\mu_p(hxC)=\mu(C)+14-s-b.
\]

τ≤μと(4)より

\[
 s+b\le\mu(C)+2.                                             \tag{5}
\]

Cが滑らかなら各直線との交点を最小にしてもs+b≥3なので不可能。
従ってCは節点三次又は尖点三次である。

各直線とCとの相異なる交点数をn_iとする。
二直線の交点がC外ならs+b=n_1+n_2+1、C上ならs+b=n_1+n_2である。
特異点を通らない直線ではn_i=1が変曲接線、2が通常接線、3が一般割線である。
特異点を通る直線では、特異接線がn_i=1、その他が2である。
節点には二つの特異接線、尖点には一つしかない。

### 2.1 節点三次の全五配置

μ(C)=1なので(5)はs+b≤3である。
両線が節点を避ける場合、節点を一線だけが通る場合、両線が通る場合の順に
n_iの上の一覧を用いると、次の五配置と一つの除外族だけが残る。

1. 二本の異なる変曲接線。
2. 節点接線と変曲接線。
3. 節点と変曲点を結ぶ割線、およびその変曲接線。
4. 二本の節点接線。
5. 変曲接線と、その変曲点を通る通常接線。

除外族は「節点接線と、節点を通る第三方向の一般直線」である。
このときΣμ=12なのでτ=12なら全局所特異点がμ=τを満たす必要がある。
しかし節点上には三つの異なる接線方向と接する二枝がある。
μ=τの平面孤立特異点は重み付き斉次であるというSaitoの判定に反する。
実際、異なる二重みの被約曲線は接線方向が高々二つであり、同じ重みなら斉次式で
全枝が相異なる直線となるからである。従ってこの全一パラメータ族は除外される。

節点三次の正規化をC*と取ると、接線の第三交点は群法でq↦q^{-2}、変曲点はμ3である。
μ3による平行移動と反転は平面自己同型へ延びるため、上の各型は一つの射影配置軌道である。
具体式で確認する場合、C=Y²Z-X²(X+Z)の正規化
[t²-1:t(t²-1):1]では、接線の第三交点のparameterは-(t²+1)/(2t)である。

### 2.2 尖点三次の全八配置

C=Y²Z-X³、尖点O=[0:0:1]、唯一の変曲点I=[0:1:0]とする。
通常点のparameter表示は[t²:t³:1]、その接線は
-3tX+2Y+t³Z=0、第三交点のparameterは-t/2である。
対角変換(t²X,t³Y,Z)が全ての非零有限parameterを動かす。

μ(C)=2なのでs+b≤4である。
両線がOを避ける場合は、変曲接線と通常接線、変曲接線とIを通る一般割線、
および二通常接線の交点が一方の接点である場合だけが残る。
二通常接線が同じ第三交点を持つことはq↦-q/2の単射性で不可能である。

一線だけがOを通る場合には、変曲接線と一般O割線、通常接線とその接点へのO割線、
通常接線とその第三交点へのO割線、尖点接線と通常接線が残る。
両線がOを通る場合は、尖点接線と一般O割線だけが残る。
したがって全部で八配置である。

ここで三つの特殊組、すなわち「尖点接線と変曲接線」、
「OIと変曲接線」、「OIと尖点接線」は、全局所式が重み付き斉次でΣτ=Σμ=13となり(4)に反する。
二本の異なる非接線O割線ではΣμ=12だが、尖点と二つの別の接線方向があるため
前節と同じμ=τ判定に反し、全補集合が除外される。
一般O割線の境界OIを黙って捨ててはいない。

## 3. 全十三配置の具体的代表と全磁気源

N_a=Y²Z-X²(X+aZ)、K=Y²Z-X³、T=-3X+2Y+Zと書く。
表の線形式をh,xとし、独立な第三線形式をyに取れば(1)の表示となる。
全二次gradient syzygyは18係数の線形系から求めた。
**全十三例でkernel次元二とSaito行列式hxCが厳密に成立する。**

| 保存ID | C | 二直線 |
|---|---|---|
| node_two_flexes | N_-3 | Z, -3X+Y+8Z |
| node_tangent_flex | N_1 | Y-X, Z |
| node_flex_secant | N_1 | X, Z |
| node_two_node_tangents | N_1 | Y-X, Y+X |
| node_flex_tangent_through_flex | N_1 | Z, X+Z |
| cusp_flex_tangent | K | Z, T |
| cusp_flex_secant_through_flex | K | Z, X-Z |
| cusp_two_tangents_incident | K | T, 12X+16Y-Z |
| cusp_flex_cusp_secant | K | Z, Y-X |
| cusp_tangent_cusp_secant | K | T, Y-X |
| cusp_tangent_cusp_third_secant | K | T, X+2Y |
| cusp_cusp_tangent_ordinary_tangent | K | Y, T |
| cusp_cusp_tangent_cusp_secant | K | Y, Y-X |

任意のR=r0h²+r1hx+r2hy+r3x²+r4xy+r5y²を持つ全源(3)を各配置で生成し、
全十のstrict Jacobi式とPfaffianをQQ上の普遍恒等式として確認した。
源を定数・線形係数だけへ制限していない。
各源の一つの固定点(r0,…,r5)=(1,2,3,5,7,11)において、全350×150 Jacobi微分と
全六源・PGL・Euler・scale方向を計算した。全十三例で実接方向の階数は48である。

入力・全syzygy・anchor・全源・全周囲証明書は
`CLASSIFICATION_TWO_LINES_PLANE_CUBIC_20261107.json`と
`CLASSIFICATION_TWO_LINES_PLANE_CUBIC_THIRD_SECANT_20261107.json`に保存した。
後者は「接点へのO割線」と異なる「第三交点へのO割線」の追加源だけを扱い、前の入力を改変していない。
実装は`code/classification_two_lines_plane_cubic.py`である。

## 4. 九つの23次元成分の支持証明

### 4.1 二次障害による六成分

次の五源では全Dの階数101、実source接方向の階数48、全kernel次元49である。
kernelから実source接空間に補な一本Rを選ぶ。
全Jacobi二次式のR²係数をDの全余核へ射影すると非零である。
表は全D列を殺す独立な有理functionalによる具体的pairingである。

| 全源 | functionalの非零項数 | pairing |
|---|---:|---:|
| cusp_flex_tangent | 70 | -2851200000 |
| cusp_two_tangents_incident | 43 | 45087997423898520000/41 |
| cusp_flex_cusp_secant | 40 | 3 |
| cusp_tangent_cusp_secant | 43 | -1240029000000 |
| cusp_cusp_tangent_ordinary_tangent | 37 | -13824 |

これは射影Jacobiの一部の式又は仮のnormalだけの計算ではない。
実48次元接方向を横切る48本の線形cutを選ぶ。
残る102次元の真の切片では、Dの101独立方向を形式陰関数で消去した後、変数は一つである。
上の余核pairingにより、その残idealは非零二次先頭項を持ち、被約支持は原点だけである。
従ってKrullの次元不等式は元の150係数空間の局所次元を48以下に抑える。
実sourceの微分階数48が与える下界と一致するので、Euler商24次元とscale一方向を引き、
PGL飽和した各既約源の閉包は全Poisson空間の**23次元既約成分**となる。

この結論は形式的なsource×normalの積分解や、正確な冪零指数を仮定しない。
全源は既約なPGL5×A6の像なので、一点の局所次元判定がその全閉包の成分性を与える。

producerは`CLASSIFICATION_TWO_LINES_PLANE_CUBIC_NORMALS_20261107.json`。
各一本のnormal、全Dを殺す疎functional、非零pairingを保存した。
`CLASSIFICATION_TWO_LINES_PLANE_CUBIC_COMPONENT_REVIEW_20261107.json`は
別のFraction実装で全五源の普遍Jacobi/Pfaffian、全六微分、D101・group47・source48・kernel49、
全Dを殺すfunctionalと非零二次pairingを再構成し、PASSした。

さらに節点型`node_tangent_flex`では、全Dの階数100、実source接方向の階数48であり、
実sourceに補な法線が二本ある。全三つの二次Jacobi係数を余核へ射影した初期空間は

\[
 \langle t_0^2,t_0t_1,t_1^2\rangle
\]

である。40,40,77項の三つの疎functionalが全Dを殺し、三二次式とのpairing行列の行列式は
-9011301580800000000000/7である。
同じ48本の線形cutと100方向の陰関数消去で残る二変数切片は原点だけを支持する。
従ってこの全六源も23次元成分を与える。
証明書は`CLASSIFICATION_TWO_LINES_PLANE_CUBIC_TWO_NORMALS_20261107.json`の同名caseである。
独立なFraction実装による
`CLASSIFICATION_TWO_LINES_PLANE_CUBIC_TWO_NORMAL_REVIEW_20261107.json`は、
全六係数族のJacobi/Pfaffian、Dの階数100、実sourceの階数48、二法線、
全三functionalによるDの消去と上記pairing行列を別に再構成し、PASSした。
この結論にも積分解、完全局所環、正確なnil指数は必要ない。

### 4.2 固定した真の切片の三次初期式による三成分

`CLASSIFICATION_TWO_LINES_PLANE_CUBIC_CUBIC_SUPPORT_20261107.json`は、
二次だけでは支持が原点にならない五つの低法線数の源について、次数三で終える有限判定を行う。
そのうち次の三源で支持原点が得られる。

| 源 | rank D | 二次初期ideal | 三次functionalの十分な値 |
|---|---:|---|---|
| cusp_tangent_cusp_third_secant | 101 | 0 | 840558307758750 t0³ |
| node_two_flexes | 100 | (t0²−16t1², t0t1−4t1²) | (t0,t1)=(4,1)で230605056000000 |
| node_flex_secant | 100 | (t0²+9t1²/4400) | qを法として−174254976000 t0t1²/11 |

全て実sourceの階数は48である。ここで用いた切片と三次式を厳密に定義する。
Nを保存した真正normal列、IをDの独立列の添字集合とする。
F+N=ker Dであり、D|span(e_i:i∈I)は単射だから

    L = span(N) + span(e_i:i∈I)

は実source接空間Fの真の補空間である。Lの消滅子が48本のliteralな線形cutを与える。
DのI列で選んだ正則minorの行集合をJとする。
切片内のD方向を陰関数で消去する際の二次補正T2(n)は、

    D[J,I] T2_I(n) = −Jac(Nn)[J]

で一意に定まる。このT2は報告の`quadratic_corrections`に各二次monomial係数として保存した。
一般にはJac(Nn)全体がim Dに入るとは主張せず、残る二次式は表のidealに含まれる。

λが全D列と全normal二次係数Qijを消すとき、λを作用させたIFT後の残式には
次数二以下の項がなく、その三次先頭項は正確にλ B(Nn,T2(n))である。
三次以降の陰関数補正はDの像なのでλによって消える。
この主張は固定した切片内の一意な補正に関するものであり、任意のkernel補正による不変性、
族の局所積分解、又はsource方向のisotrivialityを必要としない。

上表の疎functionalは順に38項、75項、13項であり、全Dと全Qijの消去も恒等的に確認した。
節点二変曲接線型の二次零点はt0=4t1の一本の直線で、三次はその上で非零である。
節点割線型ではq=t0²+9t1²/4400を法とする三次が非零定数×t0t1²であり、
q=0との共通零点は原点だけである。
尖点第三交点割線型では一変数のt0³が非零である。
従って三型とも切片の初期idealの被約支持は原点、切片そのものも零次元となる。
§4.1と同じ48cut/Krullの上界と実source48の下界により、三つの全源閉包が23次元成分となる。

独立検算に必要なfunctionalは報告の`cubic_functionals`の順に
尖点型[0]、二変曲接線型[3]、節点割線型[1]である。
最後の三次の全式は2881200 t0(3977200 t0²+2637 t1²)である。
全三次余核の再計算を要求せず、この三つと全D/源/切片の再構成で十分である。
`CLASSIFICATION_TWO_LINES_PLANE_CUBIC_CUBIC_REVIEW_20261107.json`は別Fraction実装で
全source Jacobi/Pfaffian、全D・実F・kernel、全literal102次元補切片、全TijとIFT方程式、
全保存三次係数と全functionalのD/Qij消去・pairingを再構成し、PASSした（9.05秒）。
固定切片の形式支持論も別担当が独立に受理した。

## 5. 既存成分と九成分相互の相違

本質変数数≤3は五偏微分の係数行列のrank≤3であり、五次式の空間の閉条件である。
この条件は各source閉包全体に保たれる。
二次以上の既約因子は共通の頂点P1で特異となる。
全因子が一次式へ退化した場合も、五つの超平面が同じP1を共有するためSNCではない。
従って各成分のrank4部分はSNC locusと交わらず、MPS40のいずれとも異なる。

既存の他の23次元rank4成分は、一般退化式がH+既約四次式、5H、又は三H+既約二次式である
（`ADDITIONAL_RANK4_LATE_COMPONENT_AUDIT.md`と独立三H成分の監査）。
ここでの一般退化式は相異なる二H+既約三次式なので、それらとも一致しない。
異なる次元の成分、一般rank2の六成分とは次元又は閉rank条件で区別される。
既に成分性を証明した後の比較なので、両源の稠密constructible像が同じ成分の非空開を共有するという
標準事実を用いた一般因子型の比較で十分である。

節点型の三成分は一般の三次因子が節点三次であることから、尖点型の六成分と異なる。
節点型同士では、二変曲接線、節点接線と変曲接線、節点・変曲点を結ぶ割線とその変曲接線という
三つの内在的な配置が異なる。
尖点型同士は唯一の既約三次因子、その尖点、唯一の変曲点、および二直線の接触配置から区別する。

| 全源 | 尖点を通る線 | 変曲接線 | 二直線交点のC上の位置 |
|---|---|---|---|
| cusp_flex_tangent | なし | 一つ | C外 |
| cusp_two_tangents_incident | なし | なし | 一方の通常接点 |
| cusp_flex_cusp_secant | 尖点での交点重複度2 | 一つ | C外 |
| cusp_tangent_cusp_secant | 尖点での交点重複度2 | なし | 他方の通常接点 |
| cusp_tangent_cusp_third_secant | 尖点での交点重複度2 | なし | 他方の通常接線の第三交点 |
| cusp_cusp_tangent_ordinary_tangent | 尖点での交点重複度3 | なし | C外 |

本質三変数の商平面はDを殺す定数方向から内在的に復元できるので、この表はPGL5不変である。
各sourceの一般Dは表の配置そのものであり、相異なる稠密配置は同じ成分に属し得ない。

全十三配置の被覆、fiber依存を保つSaito還元、交点budgetと全補集合の排除、相違の論理は
root担当が独立に読解・受理した。九源のKrull支持論は各独立算術担当も確認した。

## 6. 未完了の四源を保持する

残る四配置の実source次元下界は23である。全tangent上界は次の通り。

| 源 | tangent上界 |
|---|---:|
| node_two_node_tangents | 30 |
| node_flex_tangent_through_flex | 24 |
| cusp_flex_secant_through_flex | 24 |
| cusp_cusp_tangent_cusp_secant | 29 |

上界25の元の三点に対する忠実S61正規化と全134/45式の照合は、T39/T40の両graphでnegativeであった。
これは他のMPS成分への非包含や新成分性を意味しない。
上界24の残る二源の真正一法線二次・今回の固定三次は共に零であり、成分性判定には使えない。
これらでは次数を上げる反復を行わず、全sourceを保持して実変形又は配置幾何へ進む。

三次判定に先立ち、二変曲接線型の二次接錐の唯一の方向と、節点割線型の二つの共役方向に対する
literalなP+tRも全350係数で検査したが、三つともPoissonではなかった。
`CLASSIFICATION_TWO_LINES_PLANE_CUBIC_CANDIDATE_PENCILS_20261107.json`はこの限定的なnegative診断を保存する。
これは他の実変形の不存在証明ではなく、§4.2の支持証明はその診断に依存しない。
