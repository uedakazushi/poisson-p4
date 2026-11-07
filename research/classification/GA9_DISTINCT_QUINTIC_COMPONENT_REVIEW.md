# Ga9 distinct-cの新しい既約非正規quintic成分：独立検証

2026-11-07。課題E。root指定のcomponent IDは
R4_GA9_DISTINCT_IRREDUCIBLE_QUINTIC_23、
claim IDは R4_GA9_DISTINCT_IRREDUCIBLE_QUINTIC_COMPONENT_23 である。
台帳・count・Gitはこのreviewでは変更していない。

focused codeは classification_ga9_distinct_component_review.py、
証明書は CLASSIFICATION_GA9_DISTINCT_COMPONENT_REVIEW_20261107.json、
有限契約は jobs/ga9_distinct_component_review.json。
全入力は公開済みreportと既済proof noteであり、private work fileを読まない。
producer/CAS/legacy算術をimportせず、独立Fraction・疎多項式だけを使った。
約1.176秒、最大RSS19,328 KiB、CPU上限170秒、address space上限1,536 MiB。

**結論。** 表示された実曲線の非零点P1のorbit closureは、
滑らかな23次元点を一般元とする全Poisson成分である。
その一般実退化五次は絶対既約・非正規で、凍結した全prior98成分と異なる。
accepted Ga9 distinct-cの完全H+quartic源は、この成分の境界へ全て含まれる。

## 1. literal曲線、P0全源との同一性、全周囲rank

CLASSIFICATION_GA9_DISTINCT_CURVE_20261107.json のliteral cubic tensorを直接読み、
全350個の射影Jacobiがepsilonの多項式として恒等零になることを再構成した。
発見に使った六normalのIFT、有限jet、GBの完全性はこの認識の入力ではない。

\[
 t^{-1}\operatorname{diag}(1,t,1,t,t)_*P_\epsilon
       =P_{\epsilon/t}
\]

を全monomial weightで検査したので、全非零epsilon点はP1と実PGL/scalar同値である。
P0が保存anchorであることだけでなく、
元の全nilpotent source row8へaccepted canonical五parameterの零点を代入し、
全六括弧を三倍するとP0と完全に一致することも検査した。

P1の全150方向から350×150のJacobi微分Dを生成し、
全GL5作用二十五方向、Euler-vertical二十五方向、scalar一方向も別に作った。
非零102minor、非零48minorをFractionで再計算し、
Euler-vertical rank24、実group/vertical/scalar rank48、全D作用積零を確認した。
従ってDは正確にrank102であり、

\[
 \dim T_{[P_1]}\operatorname{Pois}(\mathbb P^4)
      =150-24-102-1=23,
 \qquad
 \dim\operatorname{PGL}_5[P_1]=48-24-1=23.
\]

local dimensionも23なのでP1はregular pointで、通る成分は唯一である。
そのorbit closure Cは23次元の全成分となる。
P0でのtrue cut全idealや局所形式の完全性は不要である。

## 2. QQ/Cのfactorization問題は五項modelのGauss証明で閉じる

保存された有理行列Bの可逆性と、canonical Euler代表を使った全bivectorの座標変更を
Fractionで独立に確認した。非零有理scalar倍を除く同じP1は、
座標 \((a,b,c,d,e)\) で次の九式になる。

\[
\begin{array}{lll}
 \{a,c\}=48a^2,&\{a,d\}=48ab,&\{a,e\}=32ac+b^2,\\
 \{b,c\}=-72ab,&\{b,d\}=-12b^2,&\{b,e\}=-40ad-8bc,\\
 \{c,d\}=72ad+48bc,&\{c,e\}=-48ae+6bd+32c^2,&
 \{d,e\}=72be-88cd.
\end{array}
\]

省略した \(\{a,b\}\) は零である。
全9式の変換、元の実五次のBによる変換、符号付きEuler Pfaffianを確認した結果、
このmodelの実退化式は60Fである。ここで

\[
 F=-80a^3d^2-240a^2b^2e+160a^2bcd+6ab^3d-b^4c.
\]

これはeに一次で、先頭係数は \(-240a^2b^2\)。
残りをRとすれば

\[
 R|_{a=0}=-b^4c,\qquad R|_{b=0}=-80a^3d^2.
\]

従って係数とRのgcdは \(\mathbb C[a,b,c,d]\) で一である。
Gaussの補題によりFはC上で絶対既約。
QQ因数分解結果、Galois descent、閉2+3因子locusの追加議論は使わない。

## 3. 非正規性と全reduced singular supportもnative CAS不要

\[
 F_e=-240a^2b^2,\quad
 F_c=160a^2bd-b^4,\quad
 F_d=-160a^3d+160a^2bc+6ab^3.
\]

特異点では \(a=0\) または \(b=0\)。
\(a=0\) なら \(F_c=0\) から \(b=0\)。
\(b=0,a\ne0\) なら \(F_d=0\) から \(d=0\)。
逆に次の両planeではFと全五偏微分が零であることを直接代入した。

\[
 \operatorname{Sing}(F)_{\rm red}
      =S_1\cup S_2,\qquad
 S_1=(a=b=0),\quad S_2=(b=d=0).
\]

両者はP2で、交わりはP1である。
既約三次元hypersurfaceにcodimension一の特異locusがあるのでFは非正規である。
native Singularのsingular-ideal計算はこの結論に不要。

## 4. モデルが残すrank2・nonclean rank0の具体的な補集合

以下は、全九式の確認後の直接手計算である。
S1上の接座標は \((c,d,e)\)。
Eulerをwedgeした接bivectorの係数は

\[
 c\{d,e\}-d\{c,e\}+e\{c,d\}
    =-120c^2d,
\]

従ってS1上の射影Poisson rankは一般に2である。
normal変数 \((a,b)\) に関する最低次数部は

\[
 \operatorname{in}_{I_{S_1}}F
       =80a^2d(2bc-ad),
\]

よってFのS1に沿う一般位数は3である。
これは「S1が三重divisor成分」という意味ではなく、codimension二の特異plane上の位数である。

S2上の接座標は \((a,c,e)\)。
制限したraw bivectorは

\[
 E_{S_2}\wedge(48a\partial_c+32c\partial_e)
\]

であり、射影Poisson構造としてはS2全体で零である。
raw skew matrixのrank2を射影rank2と取り違えていない。
normal変数 \((b,d)\) に関する二次symbolは

\[
 \operatorname{in}_{I_{S_2}}F
      =-80a^2(ad^2-2cbd+3eb^2).
\]

一般にはordinary doubleであるが、全symbol係数がS2上のline \(a=0\) に沿って
\(a^2\) を共通因子に持つ。
従ってS2のblowupでは、strictDとexceptional divisorの交叉に
このline上のvertical成分がmultiplicity2で現れる。
S2はまさにCLEAN rank0認識の「no vertical」仮定の補集合である。
S1はrank2かつ一般位数3なので、そこでもCLEAN仮定を誤適用しない。

## 5. 全29個の同次元MPS成分のreducibilityを実sourceから確認

diagramがdisconnectedという情報だけでは判定しなかった。
CLASSIFICATION_MPS23_DIVISOR_FACTORS_20261107.json に記録された
全29個のaccepted raw/completed generic familyを読み、
全periodと全edge parameterを残したまま、元の十括弧から
実Euler signed Pfaffianを別に構成した。
保存されたfactor積と全多項式で一致することを確認した。

| IDs | 実五次のproper factor次数 |
|---|---|
| T04–T15 | 1+1+3 |
| T18–T20 | 1+2+2 |
| T22–T31 | 1+4 |
| T32–T35 | 2+3 |

T18,T19の非monomial period分母は、可逆なperiod線形再parametrizationで
一つのLaurent変数に移して厳密に処理した。有限period標本は使っていない。
これらgeneric familyのPGL飽和閉包が該当MPS成分であることには、
既済の MPS_S61_CATALOGUE_20261103_v62.md §1–2
（MPS Theorem6.16、Corollary6.23、Example6.24による識別）を用いる。
今回の因数積恒等式だけでcomponent labelを新たに推定してはいない。

reducible quinticのconeに零を加えた集合は閉集合であり、実Pfaffian写像は多項式である。
従って全generic sourceのreducibilityから、その全component closureも同じ閉集合に入る。
境界のPfaffian零も保持される。

## 6. 凍結した全prior98との識別と全Ga9 distinct源の帰属

比較対象はimmutable CLASSIFICATION_GA9_DISTINCT_PRIOR98_20261107.json。
live台帳が後で増えても、この比較inputは変わらない。

* generic rankが異なる六成分はrankで異なる。
* rank4で次元が異なる三十六成分はfull component dimensionで異なる。
* 残る同次元23の五十六成分は、上のMPS二十九成分と非MPS二十七成分。
  非MPSのgeneric divisorも既存の証明済みfactorizationでreducible又はnonreducedである。
  両者とも全成分が閉reducible-Pfaffian locusへ入る。

Cのdense orbitの実五次は上の絶対既約Fなので、同次元の全五十六とも異なる。
これで全98とのdistinctnessが得られる。
既存二十七componentのfactorization theoremと新しく再生したMPS二十九恒等式を区別した。

最後に全非零epsilonはP1のorbitにあり、literal曲線の閉性からP0はCに含まれる。
accepted Ga9 distinct-c完全源の射影飽和閉包は次元22以下で、既約である。
P0の実group/vertical/scalar rank47から、そのorbit次元は22。
従ってaccepted全源の閉包はP0のorbit closureそのものである。
CはPGL安定な全成分なので、その全源がCへ含まれる。
この帰属にもP0での全normal IFTやhigher jetを追加仮定する必要はない。

以上は新しい全23次元成分とaccepted全sourceの帰属を支持する独立証明であり、
全Poisson分類の完成を主張するものではない。
