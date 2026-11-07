# P4内のintegral曲面のlinear normality：特異Severi定理の射程を確認した補題

2026-11-07、主義務G。以下は複素数体上のintegral・非退化な曲面
`S⊂P⁴` に適用し、smoothness、normality、ACM性を仮定しない。

**命題。** Sがlinearly normalでない場合、`(S,O_S(1))≅(P²,O_P²(2))`
であり、そのP4への埋込みは完全Veronese `v₂(P²)⊂P⁵` の同型射影である。
従って例外はsmooth projected Veroneseだけである。

## 1. 外部入力と、その仮定を実際に確認した箇所

L. Chiantini, *Lectures on the structure of projective embeddings*,
Rend. Sem. Mat. Univ. Pol. Torino 62(4) (2004), 335–388,
<https://emis.de/ft/47035> を2026-11-07に直接読んだ。
冒頭printed p.335（PDF index0）は対象を複素数体上のirreducible、reduced、
nondegenerateなprojective varietyに固定している。
**Theorem14、printed pp.369–370（PDF index34–35）** は、
1-defectiveな曲面がcone又はVeronese surfaceであると分類する。
smoothnessはglobal仮定にも当該定理にも追加されていない。
ここでは、この特異曲面を含むSeveri定理だけを外部入力とする。

独立の直接照合として、W. Niu–J. Park,
*A Castelnuovo–Mumford regularity bound for threefolds with rational singularities*,
arXiv:2203.04601v1, <https://arxiv.org/pdf/2203.04601>,
Adv. Math. 401 (2022), 108320, <https://doi.org/10.1016/j.aim.2022.108320>
の**Proposition5.3の証明冒頭、printed pp.15–16（PDF index14–15）**も読んだ。
同証明は、任意の非退化曲面S⊂P4について本命題を先に示し、
特異曲面版Severi定理としてM. Dale,
*Severi's theorem on the Veronese-surface*, J. London Math. Soc. (2)
32 (1985), 419–425を引用する。coneのvertexで接空間が全ambient空間となる点も
明記されている。Proposition5.3のstatement自体はdimension三以上の定理なので、
曲面の場合にそのstatementを誤適用してはいない。
Dale原論文はこの作業で直接精読しておらず、独立確認済み一次資料は上記二箇所である。

## 2. 射影中心はclosed secant全体を避ける

`H=O_S(1)` とし、`h⁰(H)=r+1>5` と仮定する。
Hはvery ampleなので、完全系は閉埋込み `j:S→Pʳ` を与える。
元のP4埋込みは五次元subsystem `W⊂H⁰(H)` による射影であり、
射影中心は `Λ=P^(r−5)` である。

元のWが閉埋込みを与えるため、任意の長さ二scheme `Z⊂S` に対し
`W→H⁰(Z,H|Z)` は全射である。
従って、j(Z)が張るprojective lineはΛと交わらない。
これは相異なる二点に限らず、特異点に支持される長さ二schemeにも適用する。

長さ二Hilbert scheme `Hilb²(S)` はproperである。
そのuniversal familyから得るrank二bundleとHのevaluationは各点で全射なので、
Zをそのspanへ送る `Hilb²(S)→Gr(2,r+1)` は全域で定義される。
このGrassmann像上のuniversal lineの合併はproperな像、従って閉集合である。
相異なる二点のsecantを含むので、そのclosureであるclosed secant variety
`Sec(j(S))` もこの合併に含まれる。
従って `Λ∩Sec(j(S))=∅` である。
単に相異なる二点のsecantを避けるという非閉条件に置き換えていない。

projective dimension theoremから
`dim Sec(j(S))+(r−5)<r`、従って `dim Sec(j(S))≤4`。
r≥5で非退化曲面の期待secant dimensionは五なので、j(S)は1-defectiveである。

## 3. coneの排除と唯一の例外

§1の特異Severi定理により、j(S)はcone又はVeronese surfaceである。
coneの場合、vertexの一点pからSの全ての点への直線がSに入る。
それらのtangent directionは、Sの非退化性によりambient `T_pPʳ` を張る。
従ってpでのintrinsic embedding dimensionはr≥5である。
同じ抽象scheme SはP4に閉埋込みされており、その局所embedding dimensionは四以下なので矛盾する。
これはsmoothnessによる排除ではなく、特異点のembedding dimensionによる排除である。

従ってj(S)は `v₂(P²)⊂P⁵`、r=5である。
元のS⊂P4はclosed secantを避ける一点中心からの同型射影であり、命題を得る。

## 4. rank零quintic問題への無条件の矢印

`RANK_ZERO_REGULAR_LINEARNORMAL_SURFACES.md` のEuler証明は、q=l=0なら
全rank零sectionを `I_S(2)⊗Λ²C⁵` で尽くす。
本命題により、q=0の非退化integral曲面についてl=0の追加仮定は
degree四のprojected Veronese以外では不要となる。
従って**degree五以上かつh¹(O_S)=0の全integral非退化曲面**は、
既約被約quintic退化因子を持つrank四Poisson構造のrank零曲面になれない。
normal、rational、ACM、ordinary double、又はblowupのvertical divisorに関する仮定は無い。

q>0の曲面、degree四以下の別軌道、rank二曲面、及び非正規quintic全体の分類を、
この命題だけで完了したとは主張しない。
