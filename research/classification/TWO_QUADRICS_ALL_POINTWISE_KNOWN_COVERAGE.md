# H＋相異なる二既約quadric：全点を既知有限閉和へ入れるsource合成

2026-11-08、PRIMARY E。**純粋なwhole-source合成**であり、新算術／operator／point／
minor／profile／jetを実行しない。各旧source定理、全実曲線、独立査読及び同一曲線の
有限exact readbackを区別して接続する。主稿／global台帳／Gitは変更しない。
新成分の増分は零。本稿そのものの独立aggregate査読／global採録は別の段階である。

**定理。** 複素P4上の一般階数四Poisson構造で、actual divisorが
\(D=H Q_1Q_2\)、H単純、Q1,Q2相異なる既約quadricなら、全点は以下の既知有限閉和Kへ入る。
因子rankは三・四・五、Hとの全配置、全低次／磁気係数及び全零境界を含む。
一般sourceだけ又はGの非優越をこの結論の代わりにしない。

## 1. Kの正確な既受理入力

Kは次の全閉PGL不変集合の有限和とする。最小の和又は成分数は主張しない。

* MPSの既知40成分の和。SNCに会う**全ambient成分**の尽くしという入力は
  `BETA_ZERO_COVERAGE.md` 冒頭及び `research/notes/MPS_CATALOGUE_20261103_v62_REVIEW.md`、
  named全period源は `research/notes/MPS_S61_CATALOGUE_20261103_v62.md` と
  `TWO_QUADRICS_CONFLUENCE_NAMED_MPS_PARENTS.md` に固定する。
* `BETA_ZERO_COVERAGE.md` の全Q45/零一次部分支持を被覆する既知有限閉和Kbeta。
  ここでは `legacy/v07/quadratic_classification.tex` Theorem Q22の全22源尽くしを
  **既受理分類定理として保持**し、今回再証明したとはしない。
  同noteの全22行のactual族／稠密軌道／閉性crosswalkを使い、名前の類似による同定ではない。
* CPLUS_EXT。`legacy/v08/sections/cplus_component.tex` の全成分と唯一成分点、
  `TWO_QUADRICS_PLANE_CONIC_CONFIGURATIONS.md` §5の全free-frame bundleへの包含を使う。
* 既知二成分 `R4_PLANE_TWO_CONICS_BITANGENT_23` と
  `R4_PLANE_TWO_CONICS_CONTACT4_TANGENT_23`。
  全六磁気源の成分性は `TWO_QUADRICS_PLANE_COMPONENTS.md` と二つの別作者
  `TWO_QUADRICS_BITANGENT_SEPARATE_REVIEW.md`／`TWO_QUADRICS_OTHER_TANGENT_SEPARATE_REVIEW.md`。

これらは既に定義・成分性が証明された全閉集合である。
本稿は個々の層がどの成分へ一意に属するか、交点／軌道／非被約厚みを分類しない。

## 2. 全source入口と尽くす分岐

`TWO_QUADRICS_SUPPORT.md` §§1–2の全h-Casimir liftを使い、h一で
\(\operatorname{Pf}P=\kappa q_1q_2\)、\(\kappa\ne0\)、degP≤2。
全低次項を残す逆形式の条件は
\((q_1q_2)d\alpha-d(q_1q_2)\wedge\alpha=0\)、
\(\operatorname{Pf}\alpha=kq_1q_2\)、degα≤2。
因子の互いに素性から二normal場
\(X_i=P^\#dq_i/q_i\) はaffine linearであり、
\[
 Z=\operatorname{div}P=X_1+X_2,\quad [X_1,X_2]=[X_i,P]=0,
 \quad X_1(q_2)=\lambda q_2,\quad X_2(q_1)=-\lambda q_1.
\]
両Xiは独立。第一添字／column縮約の符号は同noteに従い、
いずれを採用しても以下の零／非零・Jordan分類は同じ。

まずlambda非零／零に分ける。lambda零ではXiはh,q1,q2を全て保存し、trace零。
可換Jordan分解のjoint半単純部のspanはrank二・一・零のいずれか。
Xi二個の可換Lie代数上でJordan半単純部は線形に作用するため、これより大きいrankは無い。
この分岐はmodular場単独の有限wall表ではなく、**二つのnormal場の全支持**を分類する。
quarticの既約性を使う旧H＋既約quartic共鳴表を無条件に適用しない。

## 3. lambda非零及びjoint SS rank二：全点K

lambda非零は `TWO_QUADRICS_NONZERO_PAIR_COFACTOR_SUPPORT.md` の全matchingを使う。
全rank(4,4),(4,3),(3,3)、全6100 matchingからの三配置、
全480×90 affine線形系の完全解、全parameterのJacobi／Pfaffianを保持する。
rank四条件だけで除く一点以外の零係数も残す。
各全源のliteral全単項式SNC orbit極限にMPSのSNC全成分定理を適用して全点Kへ入る。
一般点のSNC所属だけから全点へ外挿した証明ではない。
独立入力は同NONZERO_PAIR_COFACTOR_REVIEW reportである。

lambda零・joint SS rank二では、全五weightは正確に0,±u,±v、u,v独立。
重複weightやnil partを捨てたのではなく、rank三以上の不変二次形式とtrace零から相異性を得る。
可換nil partは各一次元weight上で零となる。
全許容係数は `TWO_QUADRICS_SUPPORT.md` §§3–4のpair-zero全源
\[
 P_{xy}=C+Axy+Bzw,\quad P_{zw}=D+Exy+Fzw
\]
と四logcanonical crossで尽くされる。
全Jacobi五条件／全B=E零、B単元、E単元の補集合を保持し、
対象因子条件から必要な非零係数だけを導いた実SNC極限で全点K。
同pair-zero源のgeneric modular半単純場合もこの全源証明に含まれ、
未処理の追加resonanceを「generic」の語で捨てていない。
独立exact入力は `CLASSIFICATION_TWO_QUADRICS_PAIR_ZERO_20261107.json` と
`CLASSIFICATION_TWO_QUADRICS_MIXED_NORMAL_SUPPORT_REVIEW_20261107.json` の対応source監査。

## 4. plane全七配置：全点K

Dが本質三変数なら `TWO_QUADRICS_PLANE_WHOLE_COVERAGE.md` の全点定理を使う。
本質二変数以下ではbinary quadricが可約となるため、本対象は存在しない。
全十候補のうち三つは全域の一次syzygyでfree(2,2)と矛盾し、七つの完全Saito源が残る。
全frameGL2、全六磁気、唯一の配置parameterの全許容値／u=1/2別frameを保持する。

| whole配置 | 全被覆先 |
|---|---|
| shared_tangent_generic／contact3 | CPLUS_EXT |
| contact31_common_bitangent | 既知bitangent23成分 |
| contact4_other_tangent | 既知contact4-tangent23成分 |
| contact4_secant_at_contact | 全source actual変換でKbeta |
| contact31_chord | 真T16への全六磁気curve |
| contact31_tangent_at_simple_intersection | 真coupled T17への全六磁気curve |

旧「残二plane source未決」は作成時の履歴であり、上記全六磁気actual curve／独立receiptが
現在の陽性橋である。旧非優越や一法線診断を点wise被覆の代わりに使わない。
source入口／全配置／全Saito determinant／全frameと磁気逆は同plane aggregateの
52依存と各独立receiptに固定する。

## 5. joint SS rank一：全点K

`TWO_QUADRICS_JOINT_SS_RANK1_ACTUAL_CONFLUENCE.md` とその別作者
`TWO_QUADRICS_JOINT_SS_RANK1_CONFLUENCE_SEPARATE_REVIEW.md` を使う。
元全五型200×92必要sourceはnormal場の追加nil partを落としていない。
A1はactual一次因子、A2半単純係数非零も一次因子で対象外。
A2の半単純係数零は全nil22-image源へ戻り、§6のplane全点被覆を使う。

残B1/B2/B3の全parameterと元normalizer逆を保ったactual curveはT16又はT17。
B1の全六自由、B2のn0零／非零、B3のn1零／非零、
相異性に反するhF二乗補集合、period poleと完成項を全て保持する。
独立査読は全六括弧・実Pf・全period／frame正逆を別に手導出した。
従って旧全像≤23というG結論だけで止まっていたlaneも全点Kへ入る。

## 6. joint SS rank零：全十nil affine型を因子条件込みで戻す

Xi双方nilならZもnil。A=linearZが零なら最高次Qのdivergenceが零。
単純Hとq_i既約から \(K=\operatorname{Pf}Q=q_1|_{h=0}q_2|_{h=0}\ne0\)。
BV/Casimir恒等式 \(Q^\#dK=K\operatorname{div}Q=0\) とgeneric rank四からdK零、
非零homogeneous quarticと矛盾する。したがってA非零。
全非零nilpotent Jordan型とtranslationの零／異なるblock長topを分けた十型を使える。
この入口はK非零だけを使い、F=q1q2を既約quarticと仮定しない。
全90係数／全Jacobi／actualPfは独立NILPOTENT_SOURCES_REVIEWに固定される。

| nil型 | factor-specific whole pointwise結論 | 正確な証明 |
|---|---|---|
| 1 regular,c零 | 全T20 | 全A単元/B/C／元normalizer逆、period二次項＋nilflow＋逃げるtranslation |
| 2 regular,top4 | 空 | 全二Jacobi枝でPf零又はquadric平方、既約性を使う前の全恒等式 |
| 3 nil31,c零 | 空 | Jacobi前の全Pfにz一次因子 |
| 4 nil31,top1 | 空 | 最高quartic非零だけの全十candidate；九actual一次因子／最高部零、残一は全絶対既約quartic |
| 5 nil31,top3 | 空 | 二factor全不変二次核がx非依存、全Jacobi残源のx二乗係数非零；a零補集合もfactor対象外 |
| 6 nil22,c零 | 全plane→§4 | 全base括弧fのdiscriminant零；f非零全二閉源は重複H又は絶対既約quartic、f零は全essential三 |
| 7 nil22,top2 | 五whole源、全T16/T17/T18/T19 | 下表、元20 potential／全実normalizer逆 |
| 8 nil211,c零 | 空 | F|y零非零だけで全PDE核1,I、a0c非零による全primitive quartic；両alpha/beta零はbinary対象外 |
| 9 nil211,top1 | 全plane→§4 | 全28 potentialのF=-u c-phi+k、c非零はrank≤2因子を強制；c零はessential≤3 |
| 10 nil211,top2 | 全T17 | 元全25 potential／完成平方と全delta*a単元・全kappa（零含む）のactual曲線 |

型2/3/5/8/9/10の因子条件比較は `TWO_QUADRICS_PURE_NILPOTENT_SUPPORT.md` と
独立 `CLASSIFICATION_TWO_QUADRICS_PURE_NILPOTENT_FACTOR_REVIEW_20261107.json`。
旧表のG用次元数をEの証明とせず、空性又は全source構造だけを読む。
型4の因子型全空性と一／十の全元source逆は
`TWO_QUADRICS_PURE_NIL_1_4_10_SEPARATE_REVIEW.md`。
型6は `TWO_QUADRICS_NIL22_BASE_SUPPORT.md` §§1–5と全rank一支持独立reviewを使い、
f零の旧「未帰属二plane配置」は§4の後続全源で閉じる。
型8では旧証明中のquartic既約性から必要な仮定をF|y零非零へ正確に弱め、
二既約Q積がその仮定を満たすことを明示した因子専用還元を使う。

型七の全元20sourceは独立五族定理によって次の表を尽くす。

| whole源 | named閉親 | 全parameter／補集合 |
|---|---|---|
| rank一active-u | T17 | 全B自由、真coupled全period／inverse frame／pole相殺 |
| rank一active-v | T18 | 全c単元・全B（零含む）、period二次項でB全域 |
| rank二(eta,xi)=(0,0) | T18 | 全s（零含む）、moving親H／全h行／全projective Jacobi |
| rank二(-1/4,0) | T16 | 全s（零含む）、二parabolic frame／固定actualGL |
| rank二eta1/12,xi二乗=-1/432 | T19 | 全s（零含む）、quadric-pairを保存するactualGa²、canonical全150curve |

rank零Hessian及びrank一q=-1/2は全factor条件で空、rank二q零の二core方程式から
三coreを完全に得る。各新曲線は全spatial係数で正則、generic parameterは真completed parent、
中心は全sourceで非零。したがって閉性で全点所属を戻す。
全s単元を一orbitへ送った旧normalizerは便利な還元だが、今回の曲線はs零も直接含む。

最後のT19曲線はD102/M47のlocal診断から得た成分性ではない。
不成立だった二T19 frame及び一法線診断はそのまま履歴として残し、
別のGa²全実curveと全e/s exactreadbackによる陽性包含を使う。
型7全元20のHessian場合分けは独立 `TWO_QUADRICS_NIL22_LONG_FIVE_SEPARATE_REVIEW.md`、
新四curve及び型一／十は各WHOLE_CONFLUENCEとEXACT_READBACK及び親の独立純査読へ接続する。

## 7. 完了の論理と制限

§2のnormal場分岐は全rank四HQQを尽くす。
lambda非零は§3、lambda零はjointSSrank二→§3、一→§5、零→§6へ入る。
§6のplane像は§4で全点Kへ、その他はfactor条件から空又はnamed MPS全curve。
従って**H＋相異二既約quadricの全pointwise sourceはKに含まれる**。
ここにGの非優越だけのterminal laneは残っていない。

これは全Poisson空間の分類完了、他の退化型、H＋同一quadric二乗、複数H、
既約quadric＋cubic、非被約scheme、全軌道又は成分数上界の主張ではない。
Kへの被覆は既知成分を新たに数えない。
外部のMPS SNC theorem、既受理Q22分類、CPLUS唯一成分定理は入力として残り、
本稿の純合成や小Laurent算術がそれらを再証明したとはしない。

## 8. 保存及びhistorical bindingの区別

添付manifestは全proof／report／curve readback／source normalizerの公開依存を固定する。
旧plane manifestの四つのcurrent code/jobは互換append／binding refresh後なので、
旧manifest全hashが現pathと一致すると誤主張しない。
旧exact bytesは既公開のsource/contracts archiveから再構成してSHA一致を別記する。
旧manifest／reportを変更せず、現currentのSHAと実行時の歴史SHAを区別する。
legacy theoremは既trackedファイル又はtrackedzip内の指定memberを参照し、複製archiveを作らない。
