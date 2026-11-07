# 型7 nil22長translation第三・第四全源：相異なる二つの22次元成分の独立証明

2026年11月7日、PRIMARY E。**判定：二成分の全源包含・成分性・全prior99との相違・二者相違までPASS。**
成分数、global台帳、主TeX又はGitは本reviewでは変更していない。

公開producer入力は
`CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_QUARTIC_ORBIT_MODELS_20261107.json`、
公開七源入力は `...NIL22_LONG_QUARTIC_SOURCES_20261107.json`、
その完全有限被覆については別作者による公開 `...NIL22_LONG_QUARTIC_REVIEW_V2_20261107.json` を依存として用いる。
private work anchor/model、producer、legacy engine、CASはimportしない。
独立コード `classification_nil22_long_two_components_review.py` は、
既存の独立review用standard-library Fractionライブラリだけを再利用した。

最終証明書は `CLASSIFICATION_NIL22_LONG_TWO_COMPONENTS_REVIEW_20261107.json`、
有限契約は `jobs/hyperplane_quartic_nil22_long_two_components_review.json`。
gate PASS後の全検証は0.677秒、最大RSS20,716 KiB、CPU240秒・memory1,536 MiB上限内で終了した。
全入力hash、完全350残式、全99行と個別の識別理由を保存した。

## 1. どの全源を扱うか

全七源定理の順序を維持し、第三源はq=n14=−1、第四源はq=−1/3のrank1 V型である。
それぞれの正規parameter空間は

\[
 U_3=D(n_3-4n_6n_7)\subset\mathbb A^5_{n_0,n_1,n_3,n_6,n_7},
 \qquad
 U_4=D(n_5-n_1)\subset\mathbb A^5_{n_0,n_1,n_5,n_6,n_7}.
\]

いずれも既約である。
全V2の完全Jacobi理想両包含、逆actual core正規化、w-centeringとその全零補集合を入力とし、
正規化を戻した固定cofactor源の次元上界七、cofactor orbit十二、H四、projective scalar一から
全PGL飽和像閉包S₃,S₄はそれぞれ次元22以下となる。
この既済の**全支持上界**を一点の接空間評価で置き換えない。

新checkerは両方について、公開の全五parameter括弧を斉次化し、全350Jacobiを多項式のまま再計算した。
canonical値は第三でn3=1、第四でn5=1、他のfree parameterは零。
integer scaleはそれぞれ4、12で、全150係数が公開anchorと一致した。
parameter方向も全五本を実微分から作り、群方向への包含を確認した。

## 2. 全周囲微分・実作用と真正150次元直和

全150二次斉次方向から350×150のJacobi微分Dを構成した。
Euler外積25方向、actual GL5方向25、scalar一方向を別に生成した。
公開packetのactual GL5列はpullback微分なので、forward pushforwardの負号を明記して合わせた。
両方とも実可逆群作用の微分であり、符号を数学的な条件にしていない。

両点で厳密に

\[
 \operatorname{rank}D=102,\quad \operatorname{rank}V_{\mathrm{Euler}}=24,
 \quad\operatorname{rank}M=47,
 \quad\operatorname{rank}(M+\text{全五parameter方向})=47
\]

となる。Dの102次非零minor、全M列と全parameter列へのDの零作用、
独立な全48kernel列とその全D積をFractionで再構成した。
従ってDの上下界102が一致し、隠れたnormal方向は無い。

Mの実47列と一normal Rの結合rank48、公開48次minorは第三−191102976、第四−44236800000。
Dの独立な102標準座標列Uと合わせ、全150×150の直接行消去で

\[
 \mathbb Q^{150}=M_{47}\oplus\mathbb QR\oplus U_{102}
\]

を確認した。切片P+eR+Uuは真正な103次元線形切片である。
群方向の一部又はsourceを余分に固定して得た小切片ではない。
射影実orbit次元は47−24−1=22、射影接次元は150−24−1−102=23である。
full Jacobi schemeがsmoothだとは主張しない。

## 3. 事前有限上界を持つ完全一変数IFT

座標順(h,x,y,z,w)のC* weightsは第三(0,2,1,−1,−2)、第四(0,2,1,3,2)、
両方のtensor characterは−1、normalの相対weightは正の二である。
全Dの非零係数でrow/column weightの一致を検査し、Mの各選択列もweight homogeneousである。
選択102式のu微分minorは非零であり、その正確な値を再計算してreportに保存した。

一意な形式IFT解u(e)はC*同変なので、weight bの座標は2d=bを満たす一項ce^d以外を持てない。
uのconstant・linear termは零である。
Uの最大正weightは第三八、第四六であるから、全graphの事前次数上界は四、三。
零・負又は奇数weightの座標を、無記録の自由parameterとして残してはいない。

公開された全多項式graphを全150係数で作り直し、全102選択式へ代入すると多項式として恒等零となった。
非零IFT minorと一意性から、これは全ordersの完全形式解である。
finite jetが偶然消えたことからhigher項を外挿していない。
graphのconstant部分は正確なanchor、linear部分はR、次数二以上は全てUだけに入ることも確認した。

全350元Jacobi行をこのgraphに戻すと、第三は生成理想

\[
 (-105e^3/8,\ 4725e^4/64)=(e^3),
\]

第四は

\[
 (1485e^3/2,\ 55755e^4/8,\ -75e^2/4)=(e^2)
\]

を得る。全350の零行と全非零係数が公開packetと一致した。
各残式はweightによる一単項式であり、最小冪を持つ実際の元の行を保存したので、
理想の両包含は非零有理係数の反転だけで直ちに検証出来る。
新たなGB計算や追加nilpotence探索をしていない。

## 4. 実group×true sliceのformal productと全成分性

Euler外積の加算は、斉次二次Pに対して

\[
 E\wedge[P+E\wedge L,P+E\wedge L]=E\wedge[P,P]
\]

を満たす。GL5とscalarも全350Jacobi零schemeを保つ。
この実作用群のidentity近傍に選んだ47変数と真正切片103変数からの作用写像は、
§2の150次元直和により微分が可逆である。
従って形式inverse function theorem、又はこの点でのétale性により、
full raw Jacobi schemeの完成局所環は正確に

\[
 \mathbb C[[g_1,\ldots,g_{47},e]]/(e^3),\qquad
 \mathbb C[[g_1,\ldots,g_{47},e]]/(e^2)
\]

となる。群がschemeを保つので、引戻しはscheme-theoretically群factorと全切片の積である。
set-theoreticな被約支持だけからこのformal productを推定したのではない。

被約局所次元は47、Euler24とscalar一を除くと22。
actual22軌道閉包Cᵢは既約でこの点を通る。
この点を通る全成分は連結PGL5で保たれて軌道閉包を含み、局所次元上界22からCᵢそのものになる。
よってC₃,C₄はそれぞれuniqueな全周囲22次元成分で、dense actual orbitを持つ。

§1の全既約source像SᵢはCᵢの22軌道を含み、dim Sᵢ≤22だからSᵢ=Cᵢ。
canonical点だけでなく、全五parameterと逆正規化を戻した全H+irreducible-quartic支持が含まれる。
さらに閉包の境界も同じCᵢに入る。

## 5. 二成分の相違は既存の完全形式モデルだけで従う

もしC₃=C₄なら、二つの22次元actual PGL orbitは同じ既約成分内のopen dense orbitとなり、互いに交わる。
従って両anchorはPGL・Euler・scalarの実作用で同値になる。
しかし§4の完成局所環のHilbert–Samuel multiplicityはそれぞれ三と二であり、
実scheme automorphismに不変なので同値にならない。
同じことをgroup factorを除いた完全transverse algebraの長三・二として述べてもよい。
この識別のために新しいhigher contact、exact nil指数又は別の形式環計算を要求していない。

## 6. 全prior99との相違をimmutable catalogueから閉じる

比較は公開99commit `53da727720f053011cb53154eb789d7a6cb322be` の
`research/COMPONENTS.json`、Git blob `ad2dd611999f0ead9f9f340d0594097bb370471b` に固定した。
SHA256は `b386622e4cd64150aa515d8b901bcd88c8be462f2788197d32718c82fb241030`。
live catalogueの将来の増加に依存せず、全99行と個別理由をreview reportへ保存した。

| 既存99の群 | 個数 | 識別 |
|---|---:|---|
| generic rank二 | 6 | generic rankの違い |
| rank四・dimension22以外 | 77 | 全成分のdimensionの違い |
| dimension22・非被約D | 13 | closed非被約quintic locus |
| dimension22・1+1+3 | 1 | closed proper factor multiplication image |
| dimension22・simple H+irreducible quartic | 2 | intrinsic homogeneous cofactor Jordan型 |

新しい実Euler Pfaffianを全括弧から独立に計算すると、第三・第四はそれぞれ

\[
 D_3=2h\{(hz-wy)^2+8h^3w\},\qquad
 D_4=-6h\{(hz-wy)^2-8hw^3\}.
\]

両quarticはC[h,y,w][z]の原始二次式である。
z²係数はh²、constant termはw²y²を持つので全係数のgcdは一。
discriminantは非零scalar倍のh⁵w、h³w³で、それぞれC(h,y,w)の平方ではない。
Gaussの補題から絶対既約で、hは唯一のsimple linear factorとなる。
したがって上の13個の非被約条件、及び1+1+3の閉乗法像に入らない。
各既存成分がその閉条件に含まれることには、既存の採用済み一般因子定理を依存として用いる。

同じ22次元・H+quarticの既存二つは
`R4_NIL4_HYPERPLANE_QUARTIC_22` と `R4_NIL31_LONG_HYPERPLANE_QUARTIC_22` だけである。
この二つの公開ambient witnessの全六括弧を読み、発散を直接再計算した。
cofactorの冪rankは両方(3,2,1,0)、Jordan型(4,1)となる。
新しい二anchorの発散も直接計算し、冪rank(3,1,0)、Jordan型(3,2)、
hがdual cofactorの二乗像に入ることを確認した。

このcofactorは内在的である。
unique Hのh-Casimir代表にEuler外積を加えてなおh-Casimirなら、
hL−L(h)E=0からL=cE、従って外積自体が零となる。
よってその発散linear fieldはPGLで共役、tensor scalarで非零倍にしか変わらず、Jordan型は不変である。
既存二つもdense22 orbitを持つので、同一成分ならdense orbitが交わってこのJordan差に矛盾する。
image incidenceを全退化の閉条件としては使用していない。

これで全prior99との相違と二者の相違が閉じる。
新しい二成分のgeneric Jacobi schemeがsmoothという主張はせず、
他の五つの型7モデルのMPS帰属、全型7のE完了、全分類の尽くしもここでは主張しない。
