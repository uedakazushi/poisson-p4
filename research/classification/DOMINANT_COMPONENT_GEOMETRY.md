# 底を優越する未知成分の次元と正規退化因子

2026-11-07。未知優越成分の次元と、既存の正規退化因子分類への接続。
正規分岐と全可約・非被約因子型は処理済みであり、Gで残る一般実因子は既約・被約・非正規五次だけである。
後半の局所Lie代数補題は独立した補足で、処理済み正規分岐を再開する義務ではない。
交換子・多項式左逆・有限射影だけでは孤立残余を排除できないことは
[独立reviewの反例](NORMAL_DIRECTION_INCIDENCE_REVIEW.md)による。

## 1. 未知成分の射影次元は25、26、27のいずれか

\(W\simeq\mathbb A^{16}\)をprimitive一次部分の空間、
\(S\to W\)を元の全134式によるS61とする。以下の次元はPGL商の前である。
既知五点を分離した一般ファイバーの残余代数を\(B\)と書く。

**命題1。** \(B_{\mathrm{red}}\)が空でなければ、それに対応する全Poisson空間の
未知既約成分\(C\)が存在し、ある\(d\in\{0,1,2\}\)について

\[
 \dim C=25+d
\]

となる。ここで\(d\)は対応する一般\(W\)-ファイバーの次元である。

**証明。** 一般ファイバーの既約成分の閉包は\(W\)を優越するS61の成分を与える。
既存の[アフィン次元上界](../notes/S61_AFFINE_SURFACE_BOUND_20261102.md)により
\(d\leq2\)であり、この成分の次元は\(16+d\)である。
[C0の大域gauge監査](C0_S61_SLICE_SCALAR_AUDIT_20261106.md)の第1節を用いる。
固定点における非退化な射影評価値は5次元、primitive gaugeで除いたperspectiveは
4次元なので、全射影Poisson空間の対応成分には合計9次元が加わる。
評価値の安定化群CSp4とperspective群は連結なので、有限個の成分を置換しない。
従って\(\dim C=(16+d)+9=25+d\)である。
これはS61の55消去とBianchiの大域恒等式を入力とするC0と同じ依存を持つ。□

[既知一次像の監査](EXCEPTIONAL_BASE.md)によれば、既知の成分で\(W\)を優越するのは
T01、T39、T40のみで、その一般ファイバーが既知五点である。
従って未知\(C\)はこれらの再表示ではない。

MPSのTheorem 1.4はSNC log-symplectic locusと交わる成分を40個に尽くす。
したがって未知\(C\)はSNC log-symplectic点を一つも含まない。
一般元だけが非SNC、という弱い主張ではない。
一次資料は[MPS, arXiv:2010.08692v2](https://arxiv.org/html/2010.08692v2)、
Theorem 1.4、Example 6.24。既知40との対応は
[成分監査](KNOWN_COMPONENT_AUDIT.md)第4節を用いる。

## 2. 正規な実退化因子の分岐は既存の全点定理で閉じている

既存の `NORMAL_CLASSIFICATION` は、実際の退化因子Dが正規なら全複素閉点が
q_(5,1)、すなわちT39に属すことを証明している。
主原稿は `papers/normal_actual_divisor_recognition_ja.tex`、定理thm:classification。
第22版の定理は `legacy/v22/sections/setup.tex`、thm:normalにあり、
現在のCLAIMSでは2026-11-03の全依存scope・証明監査が完了として記録されている。

入力は正規Jacobian曲線のGorenstein認識、非活動成分の排除、E6/E8の次数制約、
四次支持の網羅性、および残る八支持の全排除である。次の既存監査を用いる。

* NORMAL_QUINTIC_RECOGNITION_AUDIT_20261103_v62
* NORMAL_ACTUAL_INACTIVE_LINE_SCOPE_AUDIT_20261103_v62
* NORMAL_ACTUAL_E6_E8_GEOMETRY_AUDIT_20261103_v62
* NORMAL_ACTUAL_RNC4_COVER_AUDIT_20261103_v62
* NORMAL_ACTUAL_EIGHT_SUPPORT_AUDIT_20261103_v62

各ノートはresearch/notes/、対応報告はreports/にある。
独立監査の341入力と、Pymの外部purely-elliptic定理を超える旧独自排除の範囲を区別する。
新しい不備が見つかったのではなく、ここでは既存結果を主分類へ接続する。

従って正規実因子の開Nについて、集合としてN⊂C_T39であり、
未知の全空間成分CはNと交わらない。もし交わればN∩CはCの非空開で、
閉集合C_T39に含まれるのでC⊂C_T39となり、成分の極大性に反する。
この包含は複素閉点での定理から被約スキームの包含へ移り、一般底への拡大後にも適用できる。
D_redだけの正規性に弱めたり、非被約Dに拡張したりしない。

## 2.1. 全次数分割を接続した現行の残余

被約可約五次の全次数分割は1+4、1+1+3、1+2+2、2+3、1+1+1+2、1^5である。
後三型は既存の全点被覆、前三型は現在の全source被覆又は射影像次元≤24で処理した。
前者の証明はHYPERPLANE_QUARTIC_NIL31_SHORT_REVIEW.md、同LONG_REVIEW.md、
NIL211_SHORT_TRANSLATION_REVIEW.md、HYPERPLANE_QUARTIC_NIL22_COTANGENT_G_BOUND.md、
HYPERPLANE_QUARTIC_NIL22_LONG_SUPPORT_REVIEW.mdと半単純/regular/最小nilの既済全支持にある。
1+1+3はTWO_HYPERPLANES_CUBIC_GENERIC_NONDOMINANCE.md、1+2+2は
TWO_QUADRICS_SUPPORT.md及びPURE_NILPOTENT_SUPPORT.mdで全型を接続する。
全非被約因子は重複H又はH+2Qで、既存の全モジュラー型非優越及び全点被覆で処理済みである。

したがって未知優越成分の一般実退化因子は、既約被約非正規quinticに限られる。
IRREDUCIBLE_NONNORMAL_ORDINARY_DOUBLE_RECOGNITION.mdは明示仮定下の一枝をT40へ認識する。
PROJECTED_VERONESE_ZERO_SURFACE_REVIEW.mdは滑らかなprojected Veronese rank0曲面を全排除する。
これらの補集合は未解決であり、G、一般五点尽くし、全分類を証明したとはしない。

## 3. 局所的な補足：正規被約退化因子とunimodular零点

\(\mathcal P\subset\mathbb P H^0(\mathbb P^4,\Lambda^2T)\)をPoisson空間とし、
\(D_\pi=V(\pi\wedge\pi)\)を退化五次因子とする。
\(\mathcal N\)を一般階数4かつ\(D_\pi\)が被約正規である部分とする。
被約Pfaffianを持つPoisson構造はlog-symplecticである。

局所体積形式\(\mu\)に対するmodular vector fieldを\(Z_\mu\)と書く。
\(\pi(p)=0\)なら\(Z_\mu(p)\)は体積形式に依らない。
実際、体積形式の変更で加わるHamiltonian vector fieldは\(p\)で零になる。
この値を\(Z_\pi(p)\)と書く。

**命題2。** 未知成分\(C\)の一般退化因子が被約正規なら、その一般元\(\pi\)は

\[
 \pi(p)=0,\qquad Z_\pi(p)=0
 \tag{1}
\]

を満たす点\(p\in\mathbb P^4\)を持つ。

**証明。** 正規log-symplectic因子の特異 locusは純余次元3のPoisson部分空間である。
4次元ではその各点を通るsymplectic leafは零次元なので、そこでは\(\pi=0\)である。
全ての特異点で\(Z_\pi\ne0\)なら、PymのDefinition 4.2とRemark 4.4により
\(\pi\)はpurely ellipticである。PymのTheorem 1.3により、この構造は
\(q_{5,1}\)、すなわちT39に属する。未知成分の一般元にはこれは起こらない。
従ってある特異点で(1)が成り立つ。□

使用する一次資料は[Pym, *Elliptic singularities on log symplectic manifolds
and Feigin–Odesskii Poisson brackets*](https://arxiv.org/pdf/1507.05668)、
Theorem 3.1、§4.1、Definition 4.2、Remark 4.4、Theorem 1.3である。
今回これらの適用仮定を本文で照合した。正規性だけでpurely ellipticとはしない。
同論文が明示する通り、\(q_{5,2}\)の退化因子は正規でないため、この分岐の対象外である。

## 4. 補足の閉incidenceと有限jet条件

\[
 \mathcal I_{\rm um}
 =\{([\pi],p)\in\mathcal P\times\mathbb P^4:
        \pi(p)=0,\ Z_\pi(p)=0\},
 \qquad
 \mathcal U=\operatorname{pr}_{\mathcal P}(\mathcal I_{\rm um})
 \tag{2}
\]

と定義する。\(\mathcal I_{\rm um}\)は閉であり、射影がproperなので
\(\mathcal U\subset\mathcal P\)も閉である。
閉性は任意のaffine座標で\(\pi(p)=0\)を課した後、一次jetの自然な縮約

\[
 T_p^*\otimes\Lambda^2T_p\longrightarrow T_p
\]

が零であるという線形条件から従う。これは局所体積形式や\(\pi\)の非零スカラー倍に
依らず、全射影空間で貼り合う。座標表示では、符号の規約を除いて
\(Z^j(p)=\sum_i\partial_i\pi^{ij}(p)\)である。

命題2と閉性から、一般退化因子が正規被約である未知成分は全体として
\(C\subset\mathcal U\)を満たす。
ただし\(\mathcal U\)は未知成分専用の集合ではなく、既知成分の境界も含み得る。
同じ証明は点ごとの包含
\(\mathcal N\setminus C_{T39}\subset\mathcal U\)も与える。

さらに(1)には以下の有限jetの必要条件がある。

**命題3。** (1)の点を原点とし、\(\pi=\pi_1+\pi_2+\cdots\)をTaylor展開する。
\(\pi_1\)は4次元unimodular Lie algebraの線形Poisson構造であり、
\(\pi_1\wedge\pi_1=0\)である。従って\(\operatorname{mult}_p D_\pi\geq3\)。
局所Pfaffianの三次斉次項を\(f_3\)とすれば
\(\pi_1^\#(df_3)=0\)である。

**証明。** Jacobiの一次斉次部分から\(\pi_1\)は線形Poissonであり、
その構造定数が\(T_p^*\)上にLie bracketを定める。
\(Z_\pi(p)=0\)は\(\operatorname{tr}(\operatorname{ad}a)=0\)を全\(a\)について
課す条件、すなわちunimodularityである。
\(f_2=\operatorname{Pf}(\pi_1)\)とする。
Poisson構造のmodular恒等式\(\pi_1^\#df_2=f_2Z_{\pi_1}\)で右辺は零である。
\(f_2\ne0\)なら\(D(f_2)\)上で\(\pi_1\)が可逆なので\(df_2=0\)。
標数零で\(f_2\)は正次数斉次式だから矛盾する。
従って\(f_2=0\)で、退化因子の次数2以下の項は消える。
全\(\pi\)の同じmodular恒等式を次数3で比較すると、
\(Z_\pi\)の定数項が零であるため\(\pi_1^\#df_3=0\)を得る。□

このLie algebra条件は必要条件である。
Lie bracketの全型を列挙しただけでは、その全高次延長や全Poisson構造を被覆したことにはならない。

## 5. 非可解一次Lie代数の分岐は全高次延長ごと排除できる

**命題4。** 滑らかな複素4次元多様体上のPoisson構造が\(p\)で零となり、
その一次Lie代数が\(\mathfrak{sl}_2\oplus\mathbb C\)なら、
\(p\)の近傍で\(\pi\wedge\pi=0\)である。
従って一般階数4のPoisson構造の零点では、一次Lie代数は必ず可解である。
この主張には退化因子の正規性やmodular vectorの消失を仮定しない。

**証明。** \(\mathfrak g=\mathfrak{sl}_2\oplus\mathfrak z\)、
\(\dim\mathfrak z=1\)とする。一次Poisson構造\(\pi_1\)に対する、係数が
\(k\)次斉次のpolyvectorの複体は

\[
 C^\bullet(\mathfrak g,\operatorname{Sym}^k\mathfrak g)
\]

というChevalley–Eilenberg複体である。
中心\(\mathfrak z\)は\(M_k=\operatorname{Sym}^k\mathfrak g\)に自明に作用するので、
複体自体が\(C^\bullet(\mathfrak{sl}_2,M_k)\otimes\Lambda^\bullet\mathfrak z^*\)
と分解する。有限次元加群\(M_k\)にWhiteheadの二補題を用いると

\[
 H^2(\mathfrak g,M_k)
 =H^2(\mathfrak{sl}_2,M_k)
  \oplus H^1(\mathfrak{sl}_2,M_k)\otimes\mathfrak z^*=0.
 \tag{3}
\]

中心が1次元で\(\Lambda^2\mathfrak z^*=0\)であることがここで必要である。

形式座標で\(\pi=\pi_1+\pi_k+\text{次数}\geq k+1\)、\(k\geq2\)まで
低次項を消したとする。Jacobiの次数\(k\)部分は
\([\pi_1,\pi_k]=0\)である。(3)より係数次数\(k\)のvector field \(X_k\)があり、
\(\pi_k=[\pi_1,X_k]\)とできる。\(X_k\)の形式flowを適切な符号で施せば
\(\pi_k\)は消え、それより低い次数を変えない。
これを全\(k\)について行う合成は、接写像が恒等である\(\mathfrak m\)-進自己同型を定める。
従って任意の全高次Poisson延長は形式的に\(\pi_1\)へ線形化される。

\(\pi_1\)は\(\mathfrak{sl}_2^*\)の3方向だけに支えられるので
\(\pi_1\wedge\pi_1=0\)。形式座標変換を戻すと、元の\(\pi\wedge\pi\)も完成局所環で零である。
滑らかな代数的または解析的局所環からその完成への写像は単射なので、実際のgermでも零となる。
この段階で形式座標変換の解析的収束は不要である。
特に既約な\(\mathbb P^4\)上では\(\pi\wedge\pi\)は全域で零となり、一般階数4に反する。

最後に、4次元の非可解複素Lie代数のLevi因子は\(\mathfrak{sl}_2\)で、radicalは1次元である。
\(\mathfrak{sl}_2\)からその1次元radicalの導分代数への作用は零なので、
Lie代数は\(\mathfrak{sl}_2\oplus\mathbb C\)である。これで非可解の場合を全て排除した。□

形式線形化は既知のLevi理論の帰結であり、新規性を主張しない。
一次資料の[Zung, *Levi decomposition of analytic Poisson structures and Lie algebroids*,
arXiv:math/0203023v2](https://arxiv.org/pdf/math/0203023)、Theorem 1.1直後のRemark 3は、
半単純代数に1次元中心を加えた場合の解析的線形化も明記している。
今回その本文を確認した。上の証明では分類に足りる形式的結論だけをWhiteheadから導いた。

この排除は任意の零点での内在的なLie代数を扱うため、S61のrank4基点を変えても成立する。
Poisson高次項を切り捨てず全次数を扱っており、一つの係数切片だけの排除ではない。
一方、可解unimodular一次部の全Poisson延長は未処理である。

## 6. 非共鳴almost-abelian一次部の局所排除

**命題5。** 一般階数4の複素代数的Poisson構造の零点pで、一次Lie代数が

    g=Ct ⋉ C³,   [t,e_i]=λ_i e_i,   [e_i,e_j]=0

と書けるとする。λ_1+λ_2+λ_3=0に加え、整数関係の格子が正確に

    {m∈Z³ : m_1λ_1+m_2λ_2+m_3λ_3=0}=Z(1,1,1)             (4)

であるなら、退化因子はpで正規ではない。
従って正規被約分岐のunimodular零点では、この一次型を全高次延長ごと排除できる。
非共鳴条件(4)をZariski開条件と呼ばず、追加整数関係を持つ型は別に残す。

**証明。** pでの形式座標を(t,x_1,x_2,x_3)とする。符号を選び、
tのHamiltonian vector field Xの線形部を

    Z=Σ_i λ_i x_i ∂_(x_i)

と書く。Jacobiより[X,π]=0である。Zは半単純である。
Xの係数次数kの各非零Z-weight項を、ad Zで割って係数次数kの座標変換で消す。
次数k未満を変えずこれを帰納すると、接写像が恒等の形式座標変換の後に

    X=Z+Y,    ord Y>=2,    [Z,Y]=0

となる。これは半単純線形部に対する形式Poincaré–Dulac操作の直接の次数別構成であり、
解析的収束を仮定しない。πも同じ座標変換で移し、[X,π]=0を保つ。

任意の非零weight μのπ成分に対して、[Y,−]はweightを保ち、係数次数を少なくとも一つ上げる。
従って(μ+ad Y)π_μ=0の最低次数の係数はμ倍となり、零以外は存在できない。
全有限jetについての帰納によりπの全項はZ-weight零である。
無限次元の固有空間の直和と直積を交換する仮定は要しない。

s=x_1x_2x_3と置く。(4)によりweight零の係数の形は厳密に

    π_(t i)=x_i f_i(t,s),
    π_(i j)=x_i x_j g_ij(t,s),    f_i,g_ij∈C[[t,s]]

となる。例えばπ_(t i)の単項式x^mではm−e_iが格子(4)に入る。
他の二座標の指数が非負なのでm−e_i=k(1,1,1)のkは非負である。
π_(i j)にも残る第三座標を用いて同じ議論を適用する。
よって形式Pfaffianは

    Pf(π)=x_1x_2x_3 F(t,s)

となる。一般階数4ならFは零でない。この形式的な退化因子の局所環は整域でない。
一方、複素有限型の正規局所環はexcellentであり、その完備化も正規である。
局所正規環は整域なので矛盾する。形式座標変換のJacobianによる単元倍は
この非正規性を変えない。□

この排除は零点での全形式延長の命題であり、固定した二次尾部だけの計算ではない。
ただし追加共鳴、非対角化可能な作用、他の可解Lie型は未処理である。
各型を列挙しただけでS61のprimitive一次像の上界が得られるとはしない。
独立担当が形式正常化の次数別構成、weight帰納、三因子の導出とnormal completionの矛盾を確認した。

**系。** 被約正規退化因子を持つP4上の一般階数4 Poisson構造の零点で、
上記の対角化可能なalmost-abelian unimodular一次型に現れる非零固有値比は、有限個に限られる。
その有限集合の値はまだ決定していない。

**証明。** 命題5によりtrace零以外に独立な整数関係が必要である。
二つの独立整数行の核は有理な一次元空間なので、固有値三組は共通の非零係数を除いてQ³に入る。
一つの固有値が非零のchartでtをrescaleし、例えばλ_1=1、λ_3=−1−λ_2とすればλ_2∈Qである。

全Poisson構造、零点p、cotangent coframe、およびλを変数とし、
一次bracketが上の型と一致する係数等式と必要な非零条件を課す。これは有限型incidenceである。
一般rank4開上の退化五次式は射影空間上の平坦なCartier因子族を与え、
被約正規fiberの条件は幾何的正規性の開条件である。
従ってこの制限も有限型であり、λ_2への像はChevalleyの定理によりconstructibleである。
Cのconstructible部分集合でQに含まれるものは有限である。
三つの非零固有値chartを合わせれば結論を得る。全λ_i=0の型は別途残る。□

この系は有限型の全incidenceに対する帰結であり、非共鳴集合をZariski開とする論法ではない。
normalityの開性・constructibility・有理比例の各段階も独立担当が確認した。
有限個の比の明示、非半単純型や他の可解型、および各全延長のprimitive一次像上界は未解決である。

## 7. Gに残る分岐

以前の草稿でG-normalとした集合は

    J_um,nor={s∈S61 : [π_s]∈(N∩U)\C_T39}

であった。第2節の既存定理によりこれは空であり、一次像の次元上界を新たに証明する必要はない。
最初の草稿がこの既存分類を接続せず未証明の主標的とした点を訂正する。
第3–6節の局所命題自体は有効だが、この分類義務を閉じるための追加計算には使わない。

未知beta優越成分はSNC点も正規な実退化因子も含まない。
従って次は、**被約だが非正規**な五次因子（可約または既約）と、**非被約**な五次因子に対し、
元134式の全支持と既知閉包を比較する。T01とT40自身がこれらの側に存在するので、
非正規性だけから未知残余を消したとはしない。
EのH+2Q全層の支持被覆は、非被約側の明示した一部分を閉じる。
一般βの全五点、他の因子型、特殊βの全被覆は未解決のままである。

## 独立確認の範囲

2026-11-07、別担当が命題1のC0次元橋、命題3のmodular恒等式と次数比較を読解確認した。
Pym一次PDFのDefinition 4.2、Remark 4.4、Theorem 1.3を別途読み、命題2の適用を確認した。
命題4についても別担当が複体の直和分解、中心1次元の必要性、全次数の形式帰納、
完成からの零性の帰還、および4次元の非可解Lie代数の帰着を読解確認した。
全134式や旧境界証明書の再計算は行っていない。
正規分岐の既存定理への接続は下記の既存独立監査を入力とし、その全計算を再実行していない。
