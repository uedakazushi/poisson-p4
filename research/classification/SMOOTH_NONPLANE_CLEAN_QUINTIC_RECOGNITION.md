# smooth nonplane clean非正規曲面：ordinary-double仮定を除いたG還元

2026-11-07、主義務G。D⊂P4をrank四Poisson tensorπの実既約被約非正規quinticとする。
S⊂SingDをsmooth integralな二次元成分とし、次の**clean条件**を仮定する。
`b:Y=Bl_S P4→P4`、`p:E=P(N_S/P4)→S`、D'をstrict transformとして、
`Γ=D'|E` の全既約divisor成分はSへdominantに写る。
codimension二のbase point上でfiber全体を含む可能性は、この条件では除いていない。

**結論。** Sがplaneでなければ、全generic transverse型を以下へ還元出来る。

* π|Sがrank零でnormal quadraticが二相異なる根を持つ場合は、Sはdegree五楕円scroll、πはT40。
* それ以外はSがsmooth projected Veroneseとなり、全πのPGL飽和像はdimension二十二以下。

従ってこの全smooth nonplane clean枝は未知G優越成分の一般元にならない。
ordinary-double仮定、normalization又は導手のsmoothnessを追加していない。

## 1. plane以外でmultiplicityは二となる

`RANK_ZERO_SURFACE_GENERIC_FIRST_JET.md` §1のelementary secant論法により、
既約quinticがnonplane曲面Sに沿いgeneric multiplicity三以上を持つことは不可能である。
従ってDのmultiplicityは正確二で、normal quadratic F2は非零である。
rank二ならgeneric transverse ordinary doubleであり、既受理のclean認識に入る。

具体的にはπ|S generic rank零なら
`IRREDUCIBLE_NONNORMAL_CLEAN_RANK_ZERO_RECOGNITION.md` のT40認識を適用する。
generic rank二なら `IRREDUCIBLE_NONNORMAL_RANK2_CLEAN_SURFACE.md` によりprojected Veroneseへ還元し、
`PROJECTED_VERONESE_TANGENT_NONDOMINANCE.md` の全source dimension二十二上界を適用する。
これらは横断nodeに関する既証明であり、以下は残るsquare normal quadraticだけを処理する。

## 2. square quadraticとclean条件はPicard群でHの二可除性を強制する

ξ=O_E(1)、H=O_S(1)とする。
DのS沿いmultiplicityは二なので

`[Γ]=2ξ+5p*H`。

normal quadraticが一般点で平方なら、Γは唯一のhorizontal prime divisor Γ0をmultiplicity二で含む。
clean条件により他のvertical divisorは無いので、**divisorとしてΓ=2Γ0**。
EはsmoothなのでΓ0はCartierであり、relative degree一から

`[Γ0]=ξ+p*B`、あるintegral line bundle B∈Pic(S)

と書ける。Picard groupのprojective bundle分解から

`2B=5H`、従って `H=2A, A=B−2H`。                     (1)

Γ0→Sが全域でsectionとなることは使わない。
isolated base point上にfiber全体があっても、上のCartier divisorとPicard classの議論は有効である。
また非零scalar residue、Gysin又はnodeの二枝を使っておらず、square normal coneのhigher transverse型を全て含む。

## 3. degree六以下のsmooth H=2Aはprojected Veroneseだけ

既存のnormalization-length degree上界からd=H²≤6。
Aはampleで、(1)から `d=4A²`。従ってd=4、A²=1である。

もしS⊂P3ならsmooth quartic hypersurfaceなのでK_S=0。
RRで `χ(A)=χ(O_S)+A²/2` がhalf-integerとなり矛盾する。
よってSはP4で非退化である。
一般H曲線のCastelnuovo上界はgenus≤1。
adjunctionとRR parityから

`g_H=3+K_S A, K_S A odd`。

従ってg_H=0、K_S A=−3。
Kodaira vanishingとrestriction sequenceによる
`H¹(O_S)→H¹(O_H)` の単射からq(S)=0。
K_S H<0なのでp_g(S)=0、従ってχ(O_S)=1。
RRとH⁰(K_S−A)=0によりh⁰(A)≥3を得る。

|A|に非零fixed divisorとmoving divisorがあれば、各々のA次数は少なくとも一なのでA²≥2となり矛盾する。
base pointがあれば解消後のfree system A'は
`A'²=1−Σm_i²≥0` より一つのsimple base pointだけを持ち、A'²=0となる。
このsystemのimageはcurveで、h⁰≥3よりprojective degree二以上。
従って一般Aは二つ以上の非零有効curveの和となり、再びA²≥2で矛盾する。
よって|A|はbase point freeである。

Aはample、A²=1なので、そのmorphismはfinite degree一でdegree一のsurfaceへ写る。
imageはP2であり、finite birational morphism to normal P2は同型である。
従って `S=P², A=O(1), H=O(2)`。
元のP4埋込みはsmooth projected Veroneseに他ならない。

## 4. Poisson rankの両方を処理する

square normal quadraticの場合の§2–3はπ|Sのrankを一切使っていない。
従ってrank零とrank二の両方でprojected Veroneseを得る。
rank零sourceにはrank四tensorが無いという既済の全source証明がある。
rank二を含む全PV-invariant sourceも、そのPGL飽和projective像がdimension≤22と独立計算・review済みである。
従ってどちらも未知25〜27次元G優越成分の一般元から排除される。

以上により、normal quadraticのrank二と一の全場合を処理した。
nonplaneではzero normal quadraticは§1により存在しない。
このためgeneric ordinary doubleでないという旧clean補集合は、smooth nonplaneでは残らない。

## 5. planeと残る正確な補集合

Sがplaneで、generic multiplicityが二かつcleanなら同じ二分岐を使える。
square quadraticでは(1)がO_P²(1)の二可除性を強制して不可能。
二相異なる根では旧rank零T40又はrank二PV認識がplaneと両立しない。
従ってこのplane subbranchにもrank四・既約被約quintic点は無い。

planeのgeneric multiplicity三以上はsecant議論の補集合であり、ここで排除しない。
smooth Sでもvertical divisorを持つ場合、及びsingular Sは引き続き残る。
特にvertical residueが高位極を持つ既知の実例と矛盾する主張はしていない。
この全clean枝のG排除は全E帰属又は既約非正規quintic全体の分類を意味しない。
