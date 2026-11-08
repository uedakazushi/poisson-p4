# nil32・rank二 CI22 の全十二源：九つの全支持と既知成分への包含

2026-11-08、PRIMARY E-FINITE-SURFACE。
**全rank四Poisson点を既知MPS又は既知H+2Q支持へ覆う。実既約被約quintic枝は全空。**
ここはrank二full12-sourceであり、先行rank零six-sourceとは座標・係数・証明を混同しない。

## 1. 正確な全sourceと九chart

`NONNORMAL_CI22_RANK_TWO_SEPARATE_REVIEW.md` §5のnil32モデルを使う。

\[
 q_1=2x_0x_2+x_1^2+2x_3x_4,\qquad
 q_2=2x_1x_2+x_4^2.
\]

元full150 tangency／intrinsic conormal kernelのEuler商がp0,…,p11の全十二源。
全十七quadricsと全350Jacobi・全126actualPfaffianは凍結された
`CLASSIFICATION_NONNORMAL_CI22_TANGENT_20261107.json`に保存され、別作者再構成済みである。
全閉点の逆caseは、先行reviewの **九つの必要map** を使う。
途中で反転したA=p10、P=p4、B=p8、f=p5、g=p11の零補集合もその九chartに戻る。
以下、各chart内のA,B,G,Z,H,Fは当該保存mapの記号であり、異なるchartの同名記号を同じ座標と見なさない。

## 2. p10=A非零の唯一compatibilityを戻す

第一mapでは A=p10、B=p8、G=p11、Z=p2、H=p7。
全p座標は先行mapの有理式で尽くされ、分母はAだけである。
このmapは五係数のnecessary supersourceであって、全部がPoissonとはしていなかった。

元17-row basisの零始まり12番は、元350係数に対し

\[
 J^{\rm basis}_{12}=-J_{41}/2+3J_{27}/4-J_{46}/2.
\]

第一mapでの値は K/(192A²)、ただし

\[
 K=384A^3H+192A^2BZ+160AB^2G+11B^4.
\]

従ってPoisson点では、既に単位であるAだけを使って正確に

\[
 H=-\frac{192A^2BZ+160AB^2G+11B^4}{384A^3}.           \tag{1}
\]

と戻る。他の非零basis residualもKの倍数だけ。
(1)を戻した全350係数は恒等零であり、残る全domainは
Spec C[A,A^{-1},B,G,Z]。新しいgeneric条件又はsource unitは加えない。

T3=x0 x4²−x1³−2x1x3x4と置くと、全actual quinticは

\[
 F=\frac{Q\,C}{3072A^4},\qquad
 Q=12A^2q_1+(8AG+B^2)q_2,                            \tag{2}
\]

\[
\begin{aligned}
C={}&384A^4T_3+192A^3B(q_1x_4-q_2x_3)
       +48A^2B^2(q_1x_2-q_2x_1)\\
 &-q_2\{(384A^3Z+256A^2BG+8AB^3)x_4
            +(192A^2BZ+128AB^2G+7B^4)x_2\}.
\end{aligned}
\]

これは保存された全126係数の恒等式である。
Cのx0 x4²係数384A⁴は非零なのでCは零多項式にならない。
q1の非退化行列Q1に対しQ1^{-1}Q2はJordan(3,2)のnilpotent。
従ってq1+λq2は全λで非退化、A非零よりQは常にrank五のsmooth quadricである。
式(2)だけでactual既約quinticは全空になる。

さらに、この**第一mapの全点はMPS閉和に入る**。
domainは既約で、Cが既約・Qと相異なる開集合は非空。
実際B=G=Z=0ではC=384A⁴T3、Q=12A²q1。
T3はx0についてprimitiveな一次式で、係数x4²と
−x1³−2x1x3x4は互いに素なので、C上既約である。
q1とは次数も異なる。
これはgeneric性を確認する手証明であり、新しい点probeではない。

その稠密開全体に既受理 `SMOOTH_QUADRIC_CUBIC_SUPPORT.md` の全smoothQ+C theoremを適用すると、
像はMPS40の有限閉和に入る。
fullsource mapはA単位domainでregularであり、同閉和の逆像も閉だから、
**Cが可約又はQと衝突する全境界も含め**第一map全体が同じ閉和に入る。
normalization-flat、familyのflatness、generic kernelの特殊化を仮定しない。

## 3. 残八mapの全actualPfaffian

既受理の九chartのラベルをそのまま用いる。
各mapの全十二座標・分母・零補集合は同名NINE_SUPPORT reportへ保存した。

| wholechart | actual signed Euler Pfaffian |
|---|---|
| p10_nonzero | §2のQ*C/(3072A⁴) |
| p10_zero_p4_nonzero | 0 |
| core_line0_nonzero | 0 |
| core_line1_nonzero | 0 |
| core_line2_nonzero、B非零 | 5x2(3Bq1−4Aq2)²/144 |
| core_zero_f_nonzero_g_minus_f | 0 |
| core_zero_f0_g_nonzero | 0 |
| all_core_zero_final_factor0 | 0 |
| all_core_zero_final_factor1 | −Z(B+Z)x2 q2² |

残八mapは元350係数が全て恒等零であり、necessary supersourceに非零compatibilityを残したまま
全点をPoissonと呼んだのではない。

core_line2ではBが元のchart unitなので3Bq1−4Aq2は必ずsmooth rank五quadric。
全rank四点は既済の H+2smoothQ 全支持に入る。
最後のmapではZ(B+Z)非零ならrank三Q=q2に対する H=(x2=0)+2Q。
Hはvertex線を含む接配置であり、全rank三H+2Q theoremの対象である。
Z(B+Z)=0の全補集合はrank二以下。
他の六mapもactualF零でrank二以下となる。
二乗quadricのsourceとgeneric Q+C sourceを混ぜず、それぞれの既受理全点定理を使う。

従って、全十二sourceのrank四Poisson点は **既知MPS40又は既知H+2rank三Q成分の閉和**へ入る。
その中にactual被約既約quinticは無い。
全mapのsource parameterが零になる境界、rankzero restriction、nonreduced R、vertical条件を捨てていない。
rank二以下は既存の全支持義務へ戻し、新しいrank四componentと数えない。

## 4. 既知全CI22分類との接続

normal CI22のrank二全空性、regular22_1のactual H+2Q factor、
singularM1の元Jacobi ideal内H factor、cone_nodeのactualD零、cone_cuspのwhole Casimir／D零と合わせ、
**integral CI22を非正規曲面支持に持つrank四・actual既約被約quinticは全空**となる。
全normal／nonnormal surface orbitは先行reviewの全symmetric-pencil分類で尽くされる。
rank零nil32六源の計算を、このrank二全十二源の代わりには使わない。
この接続はprojected Veronese、cubic scroll、他の全E又は全分類完了を意味しない。
新componentの増分は零。

## 5. 有限証明書と停止方針

`jobs/ci22_nil32_ranktwo_whole_factor.json` をgate後、最初に元17-rowspaceだけで全Pfを還元した。
そのglobal residualに十分なfactorは無く、これはpositiveとして使わず不変に保存した。
次に同じ有限契約で**既に独立被覆済みの九mapだけ**を評価し、第一mapの唯一Kを(1)で戻した。
新profile、点、minor、Jacobi operator、parameter GB又はjetは無い。
全九map・全350・全126・全零補集合のfactor readbackは
**1.189082秒、RSS53,832KiB、全exact PASS**。

global residual／全original-row liftsは
`CLASSIFICATION_CI22_NIL32_RANKTWO_WHOLE_FACTOR_20261108.json`、
全九actual sourceとfactorは
`CLASSIFICATION_CI22_NIL32_RANKTWO_NINE_SUPPORT_20261108.json`。
各実行のinline source、全map、入力SHAを保存し、旧code、旧report、currentjob、
主稿、global台帳及びGitは変更していない。
