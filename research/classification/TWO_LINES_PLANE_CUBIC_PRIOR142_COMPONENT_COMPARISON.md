# plane二H＋cubic三閉像：全prior142との比較と最後のT37義務

2026-11-08、PRIMARY E、`/root/log_manuscript`。純whole-parent監査。
新operator、point、minor、profile、局所変形jet又はCASは無し。
既known142成分、全三source及び旧prior126比較を入力にする。
本稿を保存する時点ではT37のproper bounded-leaf theoremを別作者が検証中なので、
その採録前に新成分数を加えない。

## 1. 三source閉像は既約で正確23次元

三つのactual Dは、相異二Hをh,xとし、cubicを順に

\[
 h y^2-x(x-h)^2,\qquad h y^2-(x+h)^3,\qquad h^2y-(h-x)^3
\]

とする。全R6磁気sourceのPGL5飽和閉包をX1,X2,X3と置く。
元sourceのPfは全磁気で非零、PGL5×A6は既約なので各Xiは既約閉集合である。

`TWO_LINES_PLANE_CUBIC_SUPPORT.md` の全Saito源を用いる。
固定したprojective Dの全tensorは、二次Saito場二本のGL2 frame四係数と任意R6で尽くされる。
mixed係数がframeを、fiber pairの定数base quadratic係数がRを一意に復元するので、
scalar商を取った全fixed-D fiberは正確9次元である。
反転した特定のframe entryを零へ外挿しているのではなく、det単元の全GL2を使う。

essential三変数subspaceの選択はGr(3,5)の六次元。
各marked cubic／二線配置のbase PGL3 stabilizerは零次元なので、D軌道は正確14次元。
これは全配置noteの有限軌道分類又は全logarithmic infinitesimal比較による。
PGL5によるD軌道の全fiberは上記9次元と同型だから

\[
 \dim X_i=14+9=23. \tag{1}
\]

旧全source微分F48も下界23を与えるが、(1)の証明に新しい微分やminorを要しない。
この段階では成分性又はsmoothnessは主張していない。

## 2. prior126と追加16を漏れなく合成する

旧exact prior126は
`inputs/two_lines_plane_cubic_prior126_parent_projection_20261108.json`、
全whole-parent比較の証明は `TWO_LINES_PLANE_CUBIC_KNOWN_PARENT_REDUCTION.md` に固定されている。
旧「全126」に追加された16個は、三つの23次元既約quintic成分と十三のmarked multiple-H25成分である。
添付新projectionは全142 IDを一回ずつ列挙し、旧126との集合差を明示する。
142と旧126を同じ台帳と呼ばない。

| prior142の全group | 個数 | Xi全体を含めない理由 |
|---|---:|---|
| 一般rank二成分 | 6 | closed rank条件 |
| rank四・次元23未満 | 18 | (1)との次元比較 |
| MPS23成分 | 29 | Xiの全Dはessential数≤3という閉条件、親にはessential5のSNC点がある |
| 非MPS23、generic D型が異なるもの | 46 | 等次元の包含なら等号。generic Dの既約性／素因子次数が異なる |
| 既知plane二H+cubic23の九配置 | 9 | 等号なら稠密constructible配置像が交わるが、marked cubic配置が異なる |
| nonreduced又は五Hという閉因子条件の高次元成分 | 20 | XiのDは被約で既約cubicを持つ |
| degree分割1+1+1+2又は1+2+2等を満たす24成分 | 7 | proper multiplication像にXiのDは入らない |
| T39、T40 | 2 | 第一二は同じ既存点のtangent上界24、第三は全finite lambdaのfaithful global graph陰性 |
| T21 | 1 | proper ordered-two-factor modular/contact flagの全source障害 |
| Ga9 cuspidal24 | 1 | proper marked quadratic-cone／binary discriminant障害 |
| T36、T38 | 2 | 全order二zeroのR22条件／zero-surface閉条件、各Xiの稠密開を排除 |
| T37 | 1 | 最後の義務。§5のbounded-leaf全閉定理が別途必要 |

和は142であり、未表示のknown laneは無い。
この表は**Xiというwhole閉像を一つの既知成分へ含める可能性**の比較である。
ある特殊磁気点と既知成分との交わりが空である、という強い主張ではない。

追加16のうち、三つの23成分は
`R4_PLANE_DOUBLE_LINE_CUBIC_M3_TRIPLE_4_21_23`、同`_1_12_23`、
`R4_PLANE_DOUBLE_LINE_CUBIC_M3_MOVING_SIMPLE_2_15_23`。
各generic actual Dは被約既約essential-five quinticなので、Xiと同じ23閉像ではない。
追加十三25成分は全rank一long-source成分である。
元の全sourceでh²がactual Dを割るので、properなmarked repeated-H incidenceの閉性により
**全component**がnonreduced-D条件を満たす。generic式だけの比較で境界を落としていない。
Xiのrank四点ではDが被約だから、それら十三componentへは入らない。
この追加16には、scheme smoothnessを主張しないM1/M11/B1/B3/V1も含む。
ここで使うのは成分性・次元・actual Dの閉条件であり、smoothnessではない。

## 3. 特に必要な高次元親の除外と証拠の射程

### T39／T40

第一二sourceの旧固定点の全350×150微分はD101。
Euler24及びprojective scale一を除いた全ambient tangentは24次元である。
25次元全componentの点の接次元が24になることは無いから、その同じ点はT39/T40外。
これでXi全体の包含を除くには十分であり、全lambdaについて新しいgraphや点を試す必要は無い。

第三は旧D96でtangent29なので、この議論を流用しない。
既独立 `TWO_LINES_CUSP_TANGENT_SECANT_GLOBAL_PARENT_SEPARATE_REVIEW.md` は、
I単元全R6を同じwhole QQ[lambda] faithful S61へ送り、
T39/T40の全134/45 graph両方に定数残差を得てideal(1)を証明した。
したがってXi全体はその両閉componentへ入らない。
この陰性をI零境界へ延長しない。

### T21及びGa9 cuspidal24

`TWO_LINES_PLANE_CUBIC_KNOWN_PARENT_REDUCTION.md` §§3–4は二producerのproper flagを
別作者として全係数・全factor順序で読解したproof headである。
T21の二必要条件はC|x0がh³の倍、C|h0がx³の倍となることを強制する。
三つのCは前者を満たさず、H交換でももう一方の条件に反する。
この障害は全磁気に依存しない。

Ga9 cuspidal24の全boundaryでは、可能なsimple Hのquartic制限が4pで、
そのpの0/1jetが零、quadratic coneがH²の倍、かつmarked binary discriminantに三重根が必要。
第一二では唯一4p可能なHのquadratic coneがZ(X±Z)なのでH²条件を破る。
第三では唯一可能Hのbinary discriminantがh²x²で、三重根条件を破る。
全Hの選択、a0零及びzero discriminantもproper incidenceに保持されている。
従って全Xiをその親へ含めない。

### T36／T38

`TWO_LINES_M3_WHOLE_T36_QUADRATIC_PENCIL_SEPARATE_REVIEW.md` は第一二の全R6、
J零／追加zero planesを含め全order二zero候補をvertex lineへ戻し、
全qの非零horizontal cubicがinvariant irreducible cubicを持つことからR22を除く。
第三について同じ全finite lambda／全qの証明は
`TWO_LINES_CUSP_TANGENT_SECANT_T36_T38_SEPARATE_REVIEW.md` にある。
I零を別に未決と残しても、第三whole Xiの親を除くには十分である。

T38全componentはprojective zero locusの次元が二以上というproper閉条件を持つ。
第一のlambda非零、第二のlambda非±1、第三のlambda非-1は全zero locusが一以下なのでこの条件の外。
例外lambdaの個別membershipは保持したまま、whole XiをT38へ含める可能性を除く。
旧uniform cellの陰性だけをwhole親の陰性へ昇格していない。

## 4. 三つのXiの相互相違

もしXi=Xjなら、各PGL5×A6からの稠密constructible像は同じ非空開で交わる。
第一の唯一の既約cubicはnodal、第二第三はcuspidalなので第一とは交われない。
第二はflex線とflexを通るsecantであり、第三はcusp接線とcuspを通る一般secant。
二線とcubicのmarked contactは内在的であり、その配置はprojectiveに同値でない。
従って三つの23閉像は相異なる。
特殊磁気boundaryの交差の不存在は主張しない。

## 5. 最後のT37条件と、全ambient componenthoodに必要な論理

第一の全(I,J)非零で非代数的leafとなるproofは
`TWO_LINES_NODE_FLEX_TANGENT_NONALGEBRAIC_LEAVES.md`。
第二の全(I,J)非零及び第三のI一/lambda非-1は
`TWO_LINES_CUSP_FLEX_AND_TANGENT_SECANT_NONALGEBRAIC_LEAVES.md`。
いずれもactual cubicの一般rank二点を通る代数的leaf surfaceが存在しないことを
finite-cover exact-differential residueで証明する。

もしwhole T37の全rank四点について、actual Dの各componentの一般点を通るdegree有界の
algebraic invariant surfaceがあるというproper Hilbert-incidence theoremを別に受理すれば、
各Xiの上記稠密開をT37外へ入れられる。
generic sourceだけからwhole componentへこの性質を外挿することは本稿では行わない。

**その条件付きでの成分性。** 以上によりXiを含むambient componentは既known142のいずれでもない。
既受理 `KNOWN142_SEVEN_CANDIDATE_GLOBAL_COVER_SEPARATE_REVIEW.md` の全閉被覆は、
採録外の全componentの射影次元を23以下と証明した。
従ってXiを含むcomponent CはdimC≤23、一方Xiは正確23の既約閉集合なのでC=Xi。
よって三Xiは相異なる全Poisson componentになる。

この最後の論理に新しいambient tangent、正規方向の高次延長又はsmoothnessは不要。
T37全閉定理を採録する前に本条件付き結論を新成分数へ反映しない。
本稿単独でpointwise親／元分類全体／exact総成分数を完成扱いしない。
