# 既知成分の分類への接続監査

2026-11-06。対象は複素数体上の射影化Jacobiスキーム、PGL5商の前。
このノートはAの部分成果である。rank2支持の全証明と独立査読は2026-11-07にRANK_TWO_SUPPORT.mdへ接続した。全rank4被覆・一般五点尽くし・特殊底被覆を証明しない。
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

## 4. rank4下界52の証拠と未監査範囲

40個はMPS Theorem1.4のSNC log-symplectic locusと交わる成分である。
今回一次本文の定理とExample6.24の区別の結論を確認した。
全40図式の対応は`MPS40_CROSSWALK_20261031.md`とその既存証明書を入力とする。
全MPS定理の証明を今回独立再証明したとはしない。

追加11個の全空間成分性・相違の証明入口は以下である。旧版の全証明書の再実行はしていない。
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
全局所成分性の実計算を独立に再検証した範囲とは区別する。

NEW25の成分性は、同じ実点を含む既約像のaffine次元下界26と、
全projective周囲の局所次元上界25を合わせる議論である。
全51既知rank4のうち次元25以上はT01,T39,T40,CA,CBだけであり、
これら五閉包は同じ実点を含まない。他の46成分とは次元で、六rank2成分とは
閉rank条件で区別する。`Q_POSITIVE_WEIGHT6_COMPONENT25_20261103_v63_REVIEW.md`の
この推論とcodeの依存の束ね方を確認した。minor・initial ideal・五つの閉包排除を
今回全て再実行したという主張ではない。

従って現行の算術は40+11+1=52、52+6=58である。
これは既存証明を入力にした下界であり、上界・Fの完成を意味しない。

## 5. 一次資料確認と残るA義務

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
残るAの荷重は、追加11族の全周囲成分性と
相互識別の依存を現在原稿へ接続すること、NEW25のfamily・周囲上界・五閉包排除の
独立証拠を最小集合へまとめることにある。本ノートだけでA全体をdoneにしない。
既知の誤りや反例を見つけたという報告ではなく、読解・接続と依存削減の成果である。
