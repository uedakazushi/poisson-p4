# 個別 CA/CB の形式的 bad loci と正規化 fiber の全 strata

2026-11-03。独立受理済みの `AMBIENT_ORIGIN_FULL_FORMAL_MODEL_20261103_v63_REVIEW` と `AMBIENT_ORIGIN_SEMINORMALIZATION_20261103_v63_REVIEW` からの可換環論的帰結。**各個別の完成局所環**だけを対象とし、全global bad locusやCA/CBの相対埋込みへ拡張しない。

## 模型と結論

\[
\begin{gathered}
R=\mathbb C[[g_{19},z_1,\ldots,z_6,\nu,w,r,s]],\quad
B=R[t]/(t^2-\nu),\quad F=\nu-s(s-r),\\
I=(z_1,\ldots,z_6,F,wr,ws,sr^2),\quad
J=(z_1,\ldots,z_6,F,w,r),\\
A=R\oplus It\subset B,\qquad A^{sn}=R\oplus Jt\subset B.
\end{gathered}
\]

Bは正則29次元のnormalizationである。Aのconductor商は

\[
C=R/I\simeq\mathbb C[[g_{19},w,r,s]]/(wr,ws,sr^2).
\]

以下の集合の等号は**被約な閉支持**について述べる。

\[
\boxed{\operatorname{Sing}(A)=\operatorname{Nonnormal}(A)
=\operatorname{NonCM}(A)=\operatorname{Nonseminormal}(A)
=\operatorname{Supp}(C).}
\]

この集合は滑らかな20次元の三枝

* W: r=s=0（critical枝）、
* V: w=s=0（critical枝。環の座標ではr軸）、
* S: w=r=0（self-pair枝）

の和であり、pairwise/triple intersectionは共通の滑らかなG19 `w=r=s=0` である。各枝ではz6=0、ν=s(s−r)も課す。

seminormalizationでは

\[
\boxed{\operatorname{Sing}(A^{sn})=\operatorname{Nonnormal}(A^{sn})
=\operatorname{NonCM}(A^{sn})=S,\qquad
\operatorname{Nonseminormal}(A^{sn})=\varnothing.}
\]

Sは正則20次元である。

## 全conductor点で非CMとなる証明

AのconductorはB内IB、A内でもその同じidealで、商はR/Iである。従ってconductor外ではA=Bが正則。conductor上ではnormalizationと一致しないため非正規であり、特異でもある。

任意のp⊃IをRのprimeとし、対応するAのprimeをPとする。A/pAは体κ(p)と二乗零のodd部分の和なので、pの上のAのprimeは一つである。よって

\[
A_P=R_p\oplus I_p t
\]

を有限なlocal R_p-moduleとして使える。このfiberの記述はA→Spec Rのものであり、後述の**normalization fiber**とは区別する。

R_pは正則、A_Pは同じ次元を持ち、その極大idealのradicalはpA_Pである。従ってA_PがCMなら、R_p-moduleとしてもmaximal CMになり、Auslander–Buchsbaumにより自由となる。すると直和因子I_pも自由で、そのrankは1なのでprincipalになる。しかしI_pのheightは9であり、principal idealにはなれない。矛盾である。

したがってconductorの全点で非CMである。同じ証明をJに適用すると、A^snのconductor上の全点で非CMとなる。Jはheight9のprimeで、A^snのconductor商は `R/J=C[[g19,s]]`。その外ではBと一致するため、上の三つのbad locusの等号が得られる。

## Nonseminormal locus：正確なannihilator

\[
A^{sn}/A\simeq(J/I)t=(w,r)C\,t.
\]

Cのmonomial表示から

\[
\operatorname{Ann}_C(w)=(r,s),\quad
\operatorname{Ann}_C(r)=(w,rs),\quad
\operatorname{Ann}_C(w,r)=(rs)=\operatorname{Nil}(C).
\]

例えば元の一意表示 `c(g)+w f(g,w)+r h(g,r)+s k(g,s)+rs l(g,s)` をw,rでそれぞれ掛ければ、最初の二つの等式が直接分かる。idealの交わり `(r,s)∩(w,rs)` は、ambientでは `(wr,ws,rs)` であり、C内では `(rs)` である。

従ってA^sn/Aの支持はconductorの**全被約支持**。非零の有限subintegral extensionが局所化後にも残る点でAはseminormalではない。conductor外ではAはnormalなのでseminormalである。これでNonseminormal(A)も三枝全体と一致する。A^sn自身のseminormal性は受理済みの完全square/cube証明を用いる。

## 正規化 fiber の完全な型

p⊃Iでodd部分Itをnormalizerへ送るとpBに含まれるので、対応するAの点のnormalization fiberは厳密に

\[
B\otimes_A\kappa(P)=\kappa(P)[t]/(t^2-\overline\nu).
\]

Cの被約三軸上でν=s(s−r)だから、以下を得る。

|target stratum|coneでの支持次元|B→Aのfiber|B→A^snのfiber|
|---|---:|---|---|
|conductor外|29のopen|同型、長1|同型、長1|
|W\G19、V\G19|各20|一台・双対数の長2|同型、長1|
|S\G19|20|二つの被約点、各長1|二つの被約点、各長1|
|共通G19|19|一台・双対数の長2|一台・双対数の長2|

Sではν=s²であり、s≠0なら二根±sは既に剰余体上で異なる。W,V,G19ではν=0なのでfiberはκ[t]/t²である。A^snについてはconductorがSだけで、同じ計算をJに行う。表の第二列は支持の次元であり、局所環のdimensionをそのまま表すものではない。

## 形式的横断型と中心の数値

各stratumの幾何学的一般点での**横断完成環**の数値は次の通りである。κはその幾何学的剰余体とする。

|stratum|Aのdim / edim / depth|A^snのdim / edim / depth|
|---|---|---|
|WまたはVの一般点|9 / 18 / 1|9 / 9 / 9（正則）|
|Sの一般点|9 / 18 / 1|9 / 18 / 1|
|G19の一般点|10 / 20 / 2|10 / 19 / 2|

Wではwが単位、Vではrが単位なのでI_pは正則9次元底の極大idealになり、Aの横断型は `κ[[u8,t²,u_i t,t³]]` である。

Sではsが単位で、I_pは8個のlinear parameterとr²で生成され、J_pは極大idealとなる。νが単位なのでnormalizerは二つの正則sheetへétaleに分裂する。Aの二sheetは長2のmeeting `κ[[r]]/(r²)` でglueされ、A^snでは被約点でglueされる。元のAのedimは18であり、二normalization tangentの和17とは異なる。20個のsmooth方向を戻すと、それぞれ38と37という既存の区別に一致する。

極大idealでの全cone完成環の数値は、Aが **dim29 / edim39 / depth21 / multiplicity2**、A^snが **dim29 / edim38 / depth21 / multiplicity2**。後者はJが正則列9本のidealであることからも得られる。一般点の横断数値はlocal Auslander–Buchsbaumと最小生成数から計算できる。

各環の**それぞれの**conductor上でmultiplicityは2、その外は正則で1である。実際、conductor点のlocal R_pにおいて、odd generator itの平方はνi²∈p²に入る。従ってpA_Pはm_(A_P)のreductionである。有限R_p-moduleのrank2とR_pの正則性からmultiplicityは2となる。同じ証明をJに適用する。

射影chartでは共通のsmooth因子g19がg18になり、coneのdimension・edim・depthと支持次元をそれぞれ1減らす。従って原点のprojective完成像はAが28/38/20、A^snが28/37/20、三bad branchは19次元、共通intersectionは18次元となる。

## 範囲と保存

これは個別CAと個別CBにそれぞれ同じ抽象的な**形式局所模型**を適用した結論であり、二つの元のambient埋込みを同時に同定していない。actual CA/CB intersection D、全global bad locus、収束するanalytic座標系は本稿の対象外である。

`code/ambient_origin_formal_bad_loci_20261103_v63.py` は既受理の全model/reviewのhashを結び、有限個のmonomial colon/lcmでannihilatorの等号を再現する。GBは使わない。結果は `reports/AMBIENT_ORIGIN_FORMAL_BAD_LOCI_20261103_v63.json`。
