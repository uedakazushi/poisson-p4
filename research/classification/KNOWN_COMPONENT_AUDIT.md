# 既知成分の分類への接続監査

2026-11-07。対象は複素数体上の射影化Jacobiスキーム、PGL5商の前。
現在の既知99成分について、主分類に必要な成分性・次元・相違の依存監査を完了した。
rank2支持はRANK_TWO_SUPPORT.md、追加11族とNEW25の全周囲証明は下記三監査へ接続する。
全rank4被覆・一般五点尽くし・全特殊底被覆は証明しない。
旧証明書の再実行と今回読んだ数学的議論を区別する。

## 1. 主分類で必要なrank2命題を弱める

以下、Pを全Poisson空間、P4を一般階数4の点からなる開集合、C1,...,C6を
Rat113,Rat122,Log1112,LPB,Aff5,Aff41の像閉包とする。必要な命題は

    P_{≤2} ⊂ C1∪...∪C6∪closure(P4)

で十分である。rank2 locus自体にちょうど8成分あることは、これより強い別の主張である。

この弱い命題の入力は次の三つである。

1. 因子的零点のないrank2構造は、上記六族、toricな二ベクトルの外積族、Aff32の
   閉包の和に入る。安定分岐には`R2_STABLE`の独自証明を使い、非安定分岐には
   Corrêa–Muniz v5本文Theorem4.1を使う。
2. 因子的零点を持つrank2構造はLog1112またはtoric外積族の閉包に入る。
3. toric外積族とAff32は、ともにclosure(P4)に入る。

入力2の証明。飽和foliationの次数dと零因子の次数eについてd+e=2なので、
(e,d)=(2,0),(1,1)で尽くす。次数0は線形射影、次数1はLPT Theorem6.2の
有理型(1,1,2)または線形引戻し型である。前二者の乗数付き二形式はLog1112の
特殊residueであり、Euler水平なresidue空間の二外積は分解可能である。
線形引戻しはH∂4∧v_Aで表される。既約な(H,A)族の稠密開H(e4)≠0では
H=x4とし、一般Aを対角化すれば二つの対角場の外積となる。よって全族はtoric外積族の
閉包に入る。変形はrank≤2内に保たれる。

入力3の証明。toric Lie代数の任意の交代二ベクトルはPoissonであり、分解可能な元を
一般の階数4交代元へ変形できる。Aff32はv14の普遍族Π_{β,t}を使う。
同族は全β,tでprojective Jacobiを満たし、一般β,tでは退化五次式が非零である。
一般βとPGL作用がAff32を稠密に覆うので、全Aff32閉包はclosure(P4)に入る。

したがって、もし有限個の閉族D_iがP4を覆うなら、その和はclosure(P4)も覆う。
上の包含からPはC1,...,C6,D_iで覆われる。この接続に、toric外積族とAff32の
rank2 locus内での成分性、各foliation族の独立な成分定理、全局所冪零構造は不要である。
ただし入力1の**支持の被覆**と入力2,3は省けない。

## 2. 六族の全空間成分性を支える証明

v14の実350×150整数Jacobi微分J、鉛直像V、実族の微分Aを使う。
dim im V=24、JV=JA=0は整数恒等式である。p=1000003で
rank J≥r、rank[V A]≥sを示す非零minorは、有理数体でも同じ下界を与える。
族は既約な実Poisson族であり、その微分階数と全空間接次元の間に

    s−24 ≤ dim(image closure) ≤ dim_local(P_cone) ≤ 126−r

が成り立つ。六例では両端が等しいので、閉点での局所環は正則、像閉包は全空間の成分である。
単なるrank2 locusの接空間を使った議論ではない。

| 族 | r | s | affine次元 | projective次元 |
|---|---:|---:|---:|---:|
| Rat113 | 89 | 61 | 37 | 36 |
| Rat122 | 97 | 53 | 29 | 28 |
| Log1112 | 97 | 53 | 29 | 28 |
| LPB | 86 | 64 | 40 | 39 |
| Aff5 | 103 | 47 | 23 | 22 |
| Aff41 | 103 | 47 | 23 | 22 |

入力は`legacy/v14/data/six_rank2_components.json`と、その生成コード
`legacy/v14/code/rank2_lower_bounds/verify_rank2_components.py`。
この回はコードの作用域、整数積のoverflow除外、鉛直像の商、実family微分と
上記の等号推論を読んだ。六行のminor値・全行列を新しく再計算したとはしない。
数値再現の既存記録は`audits/2026-10-30/v14_recheck/verification_summary.json`。

## 3. 同じ次元の六族を区別する補足証明

異なる次元の族は相異なるので、残る比較は二組だけである。

**Rat122とLog1112。** rank2構造πに接するprojective線形場の空間は
`ker(sl5 → H0(Λ3T), X ↦ X∧π)`である。その次元が1以上という条件は
πの係数に関する行列式の閉条件である。

Rat122の既存witness L=x0、Q1=Σxi²、Q2=Σ(i+1)xi²を考える。
葉の有理第一積分はQ1/L²,Q2/L²である。線形場Xが接すれば
L·X(Qj)−2Qj·X(L)=0となる。Lで割った剰余とL∤QjからX(L)=cL。
cEulerを引くとX(L)=0かつX(Q1)=X(Q2)=0である。
Xの行列をM、D=diag(1,2,3,4,5)とするとM+Mᵀ=0、DM+MᵀD=0。
相異なる固有値から全off-diagonalが零、最初の式からdiagonalも零となり、
projectiveな接線形場は存在しない。

一方、一般Log1112では三一次式をx0,x1,x2とできる。二次式をQとすると
X=(∂4Q)∂3−(∂3Q)∂4は三一次式とQを全て消し、非零である。
これはlogarithmic二形式の核に入るので葉方向の大域線形場である。
上の閉条件はLog1112の像閉包全体で成立するが、Rat122のwitnessでは成立しない。
よってこの二成分は相異なる。

**Aff5とAff41。** 因子的零点のない一般のaffine-Lie外積では、二つの大域線形場が
飽和接層を余次元1で自明化する。両層の反射性からT_F≅O²であり、
H0(T_F)はambientな二次元Lie部分代数そのものを復元する。
そのderived直線はnilpotent場Nの直線である。スカラー倍と共役ではJordan分割は変わらない。
分割(5)と(4,1)は異なるので、この二族の稠密な一般部分は共役でない。
仮に像閉包が同じなら、二つの稠密constructible像は非空開で交わるため矛盾する。

ここで使う「一般に因子的零点なし」は、既存witnessでも直接確認できる。
E,N,Hを列とする3×3 minorを考える。(5)ではNが一つのJordan鎖、
Hの対角成分が0,1,2,3,4である。行012,013のminorはそれぞれ
x0(2x0x2−x1²), x0(3x0x3−x1x2)で、共通因子はx0だけである。
行123のminorは2x1²x3−x1x2²−x0x2x3であり、x0では割れない。
(4,1)のwitnessはHの対角成分が0,1,2,3,7で、行012と行014のminorは
x0(2x0x2−x1²), 7x0²x4、行123は上と同じである。
こちらも全minorのgcdは1である。したがって両例でprojective外積は因子上で消えない。
因子を持つ条件は固定次数の積写像のprojective incidenceによる閉条件なので、
各既約な族の一般点でも同じ不存在が成り立つ。

以上は、六族の相異性について別のfoliation成分定理を引用する必要を減らす。
全rank2支持の被覆は第1節の依存を引き続き持つ。

## 4. rank4下界52の証拠と支持だけを使う監査

40個はMPS Theorem1.4のSNC log-symplectic locusと交わる成分である。
今回一次本文の定理とExample6.24の区別の結論を確認した。
全40図式の対応は`MPS40_CROSSWALK_20261031.md`とその既存証明書を入力とする。
全MPS定理の証明を今回独立再証明したとはしない。

追加11個の全空間成分性・相違の証明入口は以下である。旧版の全証明書の再実行はしていない。
全依存の読解と独立確認はADDITIONAL_RANK4_EARLY_COMPONENT_AUDIT.mdと
ADDITIONAL_RANK4_LATE_COMPONENT_AUDIT.mdで完了した。
特に「一般点が非SNC」だけでは40成分との区別にならず、閉包全体に課すproper incidence等を使う。

| 台帳IDのR4接頭辞を除く名前 | 次元 | 成分性・相違の入口 |
|---|---:|---|
| DOUBLE_SMOOTH_QUADRIC | 22 | v03 §new, prop:tangent, thm:quadric |
| TRIPLE_HYPERPLANE_QUADRIC | 22 | v05 triple_hyperplane.tex |
| CA, CB | 28,28 | v06 ambient_families.tex, thm:AB |
| H_QUARTIC_ESSENTIAL5 | 23 | v06 nilpotent22.tex, subsec:quarticcomponent |
| ZERO_LOG_5H | 23 | v06 zero_log_component.tex, thm:zero_log |
| MINIMAL_M | 23 | v07 minimal_thickness.tex, thm:minimal-full |
| A0 | 23 | v07 nil22_v7.tex, A0の全P4成分節 |
| DPM | 23 | v07 dplus_fullp4.tex, thm:dplus-full |
| CPLUS_EXT | 24 | v08 sections/cplus_component.tex |
| CMINUS_MAG | 23 | v08 sections/cminus_component.tex |

今回M,A0,DPMおよびC±の相違に使う条件を本文で確認した。
M/A0は唯一の特異平面の三次横断型が異なり、DPMは三特異平面を持つ。
C+は次元24、C−は次元23かつ三超平面と二次式という一般因数分解を持つ。
全局所成分性の実計算を今回再実行した範囲とは区別する。
精密な完成局所環を求めず、実際の族に横断する真の形式切片の零次元性だけで
成分性を得る補題を適用した。これにより普遍高次jetや正確なnil指数を
分類の追加前提にする必要を除いた。

NEW25の成分性は、同じ実点を含む既約像のaffine次元下界26と、
全projective周囲の局所次元上界25を合わせる議論である。
全51既知rank4のうち次元25以上はT01,T39,T40,CA,CBだけであり、
これら五閉包は同じ実点を含まない。他の46成分とは次元で、六rank2成分とは
閉rank条件で区別する。`Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW.md`の
この推論とcodeの依存の束ね方を確認した。minor・initial ideal・五つの閉包排除を
今回全て再実行したという主張ではない。
NEW25_COMPONENT_AUDIT.mdでは全依存を再点検し、局所次元の上界を八つの真正初期式、
二つの三次所属恒等式と一つの6角Sylvester行列式だけで再証明した。
この新しい有限計算をQQで実行し、独立担当が両閉枝の次元評価を確認した。

ここまでで従前の算術40+11+1=52、52+6=58を支える。
さらに RANK3_MIXED_COMPONENT_REVIEW.md の独立な全350式・全150微分と
非零一法線二次障害により、階数三二重二次曲面の22次元軌道閉包が全空間成分となる。
従前58との相違は次元・閉rank条件・二次因子のrank≤3・三乗一次因子の閉条件で示した。
さらに RANK3_TWZERO_COMPONENT.md の全実源49方向と三法線の真正二次障害により、
全Tw0源の飽和閉包が24次元全成分と分かった。五初期式の根基は原点であり、
独立Fraction実装が全Dと余核を再構成した。同次元の既知MPS8/CPLUSとは
閉非被約因子条件と被約点の存在で区別する。他次元の成分とは次元で区別できる。
したがって現在は40+11+1+1+1=54、54+6=60である。
これは既存証明を入力にした下界であり、上界・Fの完成を意味しない。

## 5. 一次資料確認と監査の射程

* [Corrêa–Muniz v5](https://arxiv.org/html/2207.12880v5)：本文4.1はstable分岐を
  algebraic integrabilityまでとする。stableの具体的Jacobian認識には独自証明が必要。
  本文の一般指数には誤植があり、n=4,k=2では同次次数和5から(1,1,3)/(1,2,2)を使う。
  §2.2は引用先のpartition≥5の仮定を拡張すると明記している。
  その拡張を無言で原論文の定理としない。主定理には第1節の弱い支持命題を使える。
* [LPT v3](https://arxiv.org/pdf/1112.3871)：Theorem6.2、Corollary6.3を確認。
  P4のcodimension2では次数1の有理型は二一次式と一二次式であり、もう一方は
  P3の大域線形場からの引戻しである。
* [MPS v2](https://arxiv.org/html/2010.08692v2)：Theorem1.4の対象はSNC locusと
  交わる40成分であり、全成分の上界を主張しない。

stable sheaf証明の全読解・独立査読とrank2支持の接続はRANK_TWO_SUPPORT.mdで完了した。
追加11族とNEW25についても、上の三つの監査およびその独立査読により
全周囲成分性・相違の必要依存を接続した。主原稿には横断切片補題とNEW25の短い上界を統合した。
従って現在の台帳に対するAは完了とする。新成分候補が得られた場合は、その族に対する
成分性・相違・次元を別途証明してから台帳へ加える。
全54rank4成分による支持被覆や全成分数の上界が得られたという意味ではない。

## 追加六成分の独立全周囲監査

最初の追加六成分により個別下界66、一般rank4下界60を得た。直前60からの六成分は次の通り。

| 全源 | 次元 | 成分性・独立確認 |
|---|---:|---|
| Jordan4 s=1/3,1/7,1/5,1/15 | 各22 | MULTIPLE_HYPERPLANE_NILPOTENT_SUPPORT.md、NILPOTENT_SOURCE_COMPONENTS_REVIEW.md。全D103/F47、全パラメータ恒等式。 |
| 従属三H+rank3Q | 24 | THREE_HYPERPLANES_QUADRIC_DEPENDENT.md §8。全D99/F49、二法線の真正全二次空間、独立Fraction再構成。 |
| Jordan(3,1)のn7非零全十源 | 25 | MODULAR_NIL31_COMPONENT.md。全D100/F50、全90係数からの15源と主開の全10源消去を独立確認。 |

四22成分は全非被約因子の閉条件でMPS、同次元の旧成分とは因子次数・重複度・二次rankで区別。
最後二成分は唯一の三重Hのcanonical cofactor像にそのHの一次式が入るかで区別する。
新24は閉な従属三一次因子条件でMPS/CPLUS、被約因子でTw0を区別。
新25は非被約因子でMPS、三乗一次因子の有無で旧positive-weight25を区別。
成分性を得た後の他次元成分との相違は次元で従い、六rank2成分との相違は階数からも従う。
文献上の新規性は未確認であり、全成分数の上界は得ていない。

## さらに四成分を加えた下界70

| 全源 | 次元 | 独立全周囲検証 |
|---|---:|---|
| nil31 short r=−1/3 | 22 | NIL31_SHORT_COMPONENT_REVIEW報告。全D103/V24/M47/F47と全源恒等式。 |
| nil31 long E,F | 各22 | NIL31_LONG_SOURCES_REVIEW.md。全D103/V24/M47/F47、全五源、canonical cofactor。 |
| 独立三H・三接線二次 m=0 | 23 | THREE_HYPERPLANES_QUADRIC_INDEPENDENT_REVIEW.md。全D102/F48、全源上界23、全syzygy核。 |

三22成分は稠密開軌道を持つ。二重H+既約三次の三成分（旧Jordan4、short、long）は、
唯一二重因子のcanonical cofactorのJordan型(4,1)/(3,2)と、その因子の双対像flagで区別する。
long Fは三重H+rank3二次なので旧rank4二次型と異なる。非被約因子の閉条件でMPSとも異なる。
独立三H新23は本質変数数≤3の閉条件でMPS、独立三HでCMINUS、他23成分とは一般因子型で区別。
これにより全70を個別同定し、一般rank4は64、rank2は6。全支持の上界は未証明である。

相異なる三H+既約Qの残三free配置は全六磁気方向を保つ実変形でMPSへ被覆した。
またsmoothQ+既約C3は全cofactor Jordan型の独立検算でMPSへ被覆した。
これら支持被覆は新しい成分数をさらに増やす主張ではない。

## 負重み法線と被約三次因子からの四成分：下界74

| 全源 | 次元・一般因子 | 独立支持証明 |
|---|---|---|
| nil22 A=−1 | 22、4H+L | NIL22_WEIGHTED_COMPONENT_REVIEW.md、全D102/G47/真103切片、非零六次障害 |
| nil31 short r=1 | 22、4H+L | NIL31_SHORT_WEIGHTED_SLICE_REVIEW報告、全六源正常化・全IFT・非零三次残式 |
| nil22 A=−1/3 | 22、2H+3L | NIL22_THIRD_WEIGHTED_SLICE_REVIEW報告、全D101/G47、真二法線、u²/uvと非零v⁴ |
| 二超平面+既約三次 p=±1 | 22、H+L+C3 | TWO_HYPERPLANES_CUBIC_COMPONENT_REVIEW報告、全五源・D103/F47・PGL同定 |

最初の三源は各一軌道で、法線は全て負の重みなので全形式陰関数は有限多項式となる。
真の被約横断支持が原点であることから全成分性を得た。精密な非被約厚さは主張しない。
二つの4H+Lは唯一四重Hのcofactor双対二乗像への所属で区別。2H+3Lは唯一三重因子の
cofactorがJordan(3,2)で、旧同因子型二成分の(4,1)と異なる。
最後の新成分は被約一般因子を持ち、他の22次元非MPS成分と異なる。MPSは全て次元23以上。
最後の証明点の群階数46と族階数47を区別し、一般軌道次元は未決定とする。
これにより全74成分を個別同定し、一般rank4は68、rank2は6となる。全支持の上界は未証明。

## 平面曲線十二成分と下界86

DEGENERATE_QUADRIC_Q3_AZERO_SUPPORT.mdの六conic+cubic源と、
TWO_LINES_PLANE_CUBIC_SUPPORT.mdの指定六two-lines+cubic源は、全て23次元全空間成分である。
全D100/101、実方向48、真正一二法線の二次支持と必要な三次余核を独立Fractionで再構成した。
三次判定は全kernel補正の曖昧さを消したものだけを使い、特定liftの失敗からは推論しない。
conic+cubic六源は全磁気shearのrank五商とscalar作用で稠密軌道を持つ。
two-lines+cubic六源は群方向47と族方向48を区別し、一般軌道次元は未決定のまま。

本質三変数条件は全閉包に保たれ、全MPS40のSNC点を含められないため全MPSと異なる。
同次元の従前非MPSは一般因子型1+4、五一次因子、1+1+1+2で区別し、新十二同士も2+3/1+1+3と
平面配置のnode/cusp、特異点通過、交点重複度、接線/変曲接線で区別する。
独立報告はQ3_FREE_CONIC_CUBIC_COMPONENT_REVIEW、Q3_FREE_CONIC_CUBIC_THIRD_REVIEW、
TWO_LINES_PLANE_CUBIC_COMPONENT_REVIEW、TWO_LINES_PLANE_CUBIC_TWO_NORMAL_REVIEW（全て20261107）。
従って全86を個別同定し、一般rank4は80、rank2は6。全分類上界は未証明である。

## 固定切片の三次初期式による三追加成分と下界89

TWO_LINES_PLANE_CUBIC_CUBIC_REVIEW.mdの三全源を追加した。
全D100/101と実源48、真正102次元線形切片、選択minorにより一意な陰関数二次補正を独立再構成した。
全Dと全normal二次式を消す余核により、真の三次初期式を得る。
二次初期式との共通零点が原点だけであることから、Krull不等式で全成分性を証明した。
任意kernel補正の不変性は仮定せず、実際に固定した切片の一意補正を用いる。
node二変曲接線、node割線・変曲接線、cusp第三交点割線という配置は既登録六型と相異なる。
本質三変数閉条件と一般因子型を合わせ、既存86成分との相違も確認した。
独立報告はCLASSIFICATION_TWO_LINES_PLANE_CUBIC_CUBIC_REVIEW_20261107.json。
全89成分、一般rank4は83、rank2は6。一般軌道次元と全分類上界は未証明である。


## regular nilpotentの被約1+4型22次元成分と現行下界90

HYPERPLANE_QUARTIC_NIL_REGULAR_COMPONENT_REVIEW.mdは、全150係数・350式の
D102、実群・vertical・scalar47、完全kernel48と真正一本normalの二次補正を独立再構成した。
全kernelによる補正曖昧さも微分像内に入り、三次余核9108855/2は非零である。
真正切片の被約支持が原点なので、canonical PGL軌道閉包は22次元全空間成分となる。
HYPERPLANE_QUARTIC_NIL4_COMPONENT_DISTINCTNESS.mdと独立reportは、
基準89の全次元/rank/閉factorizationを監査し、この被約1+4型が全旧成分と相異なることを示す。
実退化式はh[(2wy−z²)²−48hw³]で、四次因子は原始h一次式として絶対既約。
全第一normalized源の既約像も、canonical点での唯一成分性によりこの成分へ入る。
この一成分を基準89へ加えると下界90、一般rank四84となる。さらに次節の二成分を加える。上界と全分類は未証明。


## 二平面conicの二つの23次元成分と現行下界92

TWO_QUADRICS_PLANE_COMPONENTS.mdおよび別作者BITANGENT/OTHER_TANGENTレビューは、
全六磁気源・全150係数/350式・D101/F48/kernel49と真正48cutを再構成した。
異点共通接線型の二次余核−120295458720、四点接触・別点接線型の真の固定二次補正後の
三次余核−533355840000は非零であり、被約横断支持は原点。実source次元下界と一致して
全飽和像閉包が23次元全成分となる。局所唯一性や一般軌道次元は追加主張しない。
essential三変数の閉条件と一般因子1+2+2、相互のconic交点3+1対4により旧90と区別した。
従って個別下界は92、一般rank四86、rank二6。全上界と分類は未証明。


## nil31 image内half/thirdの二つの25成分と現行下界94

NIL31_IMAGE_ZERO_COMPONENTS_REVIEW.mdで全十源のPGL/Euler橋、全Jacobi/Pf、
全150/350行列とactual F50を別作者Fractionで再構成した。
halfはD100の接上界25と源下界25が一致する。thirdはD98・全kernel52と真正100次元cutで、
二余核初期式s²−At²,st+Bt²のB²−A≠0から根基原点を得て上界を証明した。
両全飽和像閉包は全空間25次元成分である。thirdの滑らかさや局所唯一性は主張しない。
全H³L²閉条件で既存92の同次元五成分と区別し、一意な三重超平面のcanonical cofactor
双対像flagで相互を区別した。全既存92の比較入力はPRIOR92固定証明書に保存した。
従って個別下界94、一般rank四88、rank二6。全上界と分類は未証明。


## nil22 image内二節点三次23成分と現行下界96

NIL22_IMAGE_COMPONENTS_REVIEW.mdの全6磁気源・全150/350再構成、G=F48、D101/100と
真正102次元切片を確認した。法線初期式a²、またはa²−3b²/124,ab+5b²/31は最大根基を持ち、
Krull上界が実源下界と一致する。全源閉包は各23次元全成分で、実PGL軌道23は稠密開。
全閉H²C3条件と旧5H閉条件で同次元旧成分を区別し、唯一二重Hと唯一cubic singular planeの
incidenceで相互を区別した。全旧94の比較は固定PRIOR94証明書に依存する。
従って個別下界96、一般rank四90、rank二6。全Jacobi schemeの滑らかさ・局所唯一性、
全成分上界・全分類は未証明。


## nil31長translation第一源の22成分と現行下界97

QUARTIC_NIL31_LONG_COMPONENT_REVIEW.mdは全第一五源の実affine一軌道化、全90/110源、
全150/350微分D103・Euler24・G=F47を別作者Fractionで確認した。
全Jacobi schemeのcanonical点は正則唯一22成分、実軌道22が稠密開である。
全旧96の同次元rank4十五件を、十三非被約・1+1+3・旧nil4のuniqueH cofactor flagへ尽くして
区別した。h∈im(N*)³とh∉im(N*)は同じ1+4型の二denseorbitを区別する。
第二五源も別作者のfull134/45T40 graphとD100により全T40へ被覆し、型5入口を完了した。
個別下界97、一般rank四91、rank二6。全上界・全分類は未証明。


## Ga9 c=1全十源の23成分と現行下界98

GA9_CONSTANT_ANCHOR_REVIEW.mdは全十源、全350式、D98/F48、真正102cutと全150直和を
Fraction別算術で確認した。正法線weight(1,2,2,3)から全IFT10・全Jacobi16の有限性が先験的に従い、
全78補正と179残式を再構成した。原179式からの四liftにより根基原点、Krullにより全十源閉包は23成分。
閉triple-plane条件、proper因子積像、旧四1+4成分の直接実cofactor計算で全97と区別した。
原凍結97projectionは公開済みd27のGit blobと完全一致し、metadataHEADだけがcommit前の0f0cである。
個別下界98、一般rank四92、rank二6。滑らかさ・局所一意性・全分類は未証明。


## 既約非正規quinticの滑らか23軌道成分と現行下界99

GA9_DISTINCT_COMPONENT_REVIEWは全350literal三次曲線、P1全D102/G48、全可逆座標変更と九括弧、
e一次primitive五項Fの絶対既約性・二特異P2をFractionで独立確認した。
接次元=実軌道23なので軌道閉包自身が滑らかな一般点を持つ全成分、商次元0。
全Ga9二根源はcanonical稠密22軌道の実曲線極限からその境界へ入る。
同次元MPS29の全一般period/edge源の実D因子を直接恒等式で再構成し、
非MPS27の既存閉可約因子条件と合わせ全旧98と区別した。
個別下界99、rank四93、rank二6。全上界・全E/G尽くしは未証明。


## 型7の三新成分と現行下界102

全七源V2の独立被覆定理、NIL22_LONG_TWO_COMPONENTS_REVIEWとNIL22_LONG_CURVES_REVIEWから、
第三・第四の全五源閉包は相異なるdense22軌道成分、第五は新smooth23軌道成分へ入る。
前二の全350完全正weight IFT/e³・e²と実group×true slice étale積から成分性と相互相違、
unique H cofactor Jordan32と旧二22の41から旧成分と相違する。全prior99をimmutable blobから比較済み。
第五のprimitive w-linear quinticは絶対既約で、unique stabilizer semisimple multiplicities221は旧既約五次23の11111と異なる。
新三者の相違も前二対第五はdimension差で閉じる。一般rank四96、rank二6、個別下界102。上界・全尽くしは未証明。


## 平面A3+A2 quartic成分と現行下界103

GA9_CUSPIDAL_QUARTIC_SEPARATE_REVIEWは全七源/全六磁気移送、D100/F49/真正49cutの非零二次101minorとKrull、全prior102相違を独立確認。一般orbit23、商1。全Ga9 power1/power3閉包を被覆する。下界103/rank四97/rank二6、全上界未証明。


## 平面(4,5)既約quintic成分と現行個別下界104

ESSENTIAL_THREE_QUINTIC_FIRST_SEPARATE_REVIEWは全28x18gradient基底、全350/D101/M48/真正48cutの非零二次λQ、immutableprior103との全相違を独立確認。全六mag源はactual23orbit、商dimension0。下界104/rank四98/rank二6、全上界は未証明。
