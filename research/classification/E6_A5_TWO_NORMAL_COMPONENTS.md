# E6 smooth-flex／A5 smooth-hyperflexの全二法線支持：二つの23成分候補

2026-11-07、PRIMARY E。producerの全支持証明、独立算術／全prior124監査待ち。
新しい二成分をまだglobal計数しない。

## 同じ既canonical点の完全map

両点は既MU6/MU5全sourceのcanonical representativeそのものであり、primitive integer scaleも旧値のままである。
既fullD100の非零100minorを新しい全350×150微分から再計算し、旧det値に完全一致することを確認した。
実Euler/GL/scale方向の全51列はrank48。

その48pivot行で真102次元線形cutを取り、全D kernelからF48を引いて二法線N1,N2を作った。
全48+2+100方向はraw150の直和で、そのunit determinantを保存した。
D*F=D*Ni=0、100implicit方向のD像blockは既受理のunitである。
これは小さいsliceの推測を全空間へ拡張する操作ではない。

`E6_A5_TRUE_TWO_NORMAL_OBSTRUCTION.md` のformalgroup-slice／IFTの橋に従い、
全350Jacobiの二次三列(u²,uv,v²)を取り、全D imageを消した。
全350の余核行と、全150D列を消すfunctional、全行と選択basis formsの両方向係数恒等式を保持した。
選んだ二次式だけで全mapを代替していない。

## E6：全quadratic spanが三次元

E6 smooth-flexでは、全余核quadratic spanはrank三。
三独立basis formsからu²とv²をconstant線形結合で表す二つのexact identityを保存した。
従って全quadratic idealの根基は(u,v)である。
formal true cutの被約支持は原点、全raw局所dimensionは48以下となる。
実orbit下界48と一致するので、この実23軌道閉包が全ambient componentとなる。

## A5：二独立quadraticsはP1で共通零を持たない

A5 smooth-hyperflexの全quadratic spanは次の二式で尽くされる。

\[
g_1=\frac52uv,\qquad
g_2=\frac{54}{25}u^2-\frac{186624}{7225}v^2.
\]

全homogeneous degree三Macaulay行列のdetは`-2519424/7225`である。
さらに明示的に

\[
u^3=\frac{6912}{1445}v g_1+\frac{25}{54}u g_2,
\]
\[
v^3=\frac{289}{8640}u g_1-\frac{7225}{186624}v g_2.
\]

両identityを全係数で検証した。無限遠を除くaffine root検査ではなく、全P1に共通零が無い証明である。
従ってA5でも全true cutは被約零次元、実23軌道閉包が全ambient componentになる。

両者ともtangentはraw50、実dimensionは48であり、scheme smoothnessや精密nil厚みを主張しない。
高次Jacobi jetは使っていない。

## 全六磁気と境界

全source normalizationは既MU6/MU5別作者reviewで受理済みである。
各固定marked pairについて全六磁気sourceの閉包は、今回使った一23orbitの閉包に等しい。
これは非零magnetic invariantの主開だけの主張ではなく、零invariantの全境界をregular sourceの閉性で戻している。
従って上の全componentは元の完全六磁気sourceとそのPGL飽和閉包を覆う。

全二次mapの一括計算は10.023秒、最大RSS78,380KiBで終了した。
report `CLASSIFICATION_PLANE_HYPERPLANE_QUARTIC_TWO_NORMAL_OBSTRUCTION_20261107.json` はCOMPLETE。
exact execution bytesもreportに保存し、旧17 FunctionDef ASTを変更していない。

## 全prior124及び相互の相違

priorは公開commit `8b46f29ddbb9f539abc7a026e492ee25bf8ad699` の124行全catalogueをimmutable projectionへ固定する。
既成分のgeneric rank又は全component dimensionが異なれば同じcomponentではない。
同じ23のMPS29はaccepted actual toric-SNC点でDが非零κ*hxyzw、essential rank五。
新orbit閉包全体のessential≤3というclosed条件と矛盾する。

旧same23/reducedH+quartic五源は実witnessで扱う。
ESSENTIAL5とGa9constantのactualD essential rank五は既MU6独立reviewの全gradient minorで確認済み。
MINIMAL_Mの実cofactor rank一は新rank二と異なる。
A0はD5、DPMは3A2であり、今回のquartic E6/A5とはintrinsic singularity profileが違う。
他の旧成分はclosed factor multiplication、nonreduced divisor、又は旧irreducible quintic actual点で区別する。

既採録の新μ6九源とμ5四源についても全marked profilesを保持する。
E6 smooth-bitangent旧成分だけは同じE6 quarticを持つが、H上contact partitionが(2,2)である。
今回のE6 smooth-flexは(3,1)なので、両dense23orbitsは交わらない。
A5 quarticは他十三源のA6/A4+A2/3A2/E6/A4+A1/D5/2A2+A1のいずれとも異なる。
今回二つ相互もE6 unibranch cuspとA5 two smooth tangent branchesなので相違する。

同じ23componentなら両actual23orbitsはそのdense open orbitであり、交わらなければならない。
従ってこの比較は、境界で不変か未確認のincidenceを無理にclosedと仮定するものではない。
全124行のdispositionは別auditに保存する。正式独立reviewとroot読解を経るまで計数しない。

これで対象二源のE義務は、証明が受理されればそれぞれこの全23componentへの帰属で閉じる。
Ga9 linear-w、stable二源、nilshort r=2/3その他の未完了源はこの結果だけでは処理しない。
