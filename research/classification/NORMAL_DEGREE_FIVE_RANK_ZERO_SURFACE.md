# normalな次数五rank零曲面は全てT40を認識する

2026-11-07。主義務Gの既約被約非正規quinticに対する部分定理。
本文はrootによる幾何証明であり、NORMAL_DEGREE_FIVE_RANK_ZERO_SEPARATE_REVIEW.mdで別作者が一次資料と全適用範囲を独立確認してPASSした。

## 命題と射程

S⊂P4をintegral・非退化・normalな次数五曲面とし、rank四Poisson構造πがS全体で零、
実退化因子Dが既約被約五次とする。このときSは滑らかな次数五楕円scrollであり、
πは同じ埋込みのq52のscalar倍、従ってT40に属す。
Sのsmoothness、Dの横断node性、blowup交叉のvertical成分不存在を仮定しない。
次数六、非normal S、rank二restrictionは本命題に含めない。

## 1. 一次資料の適用範囲

A. L. Tironi, *Normal projective varieties of degree 5*, arXiv:1201.4508v1,
<https://arxiv.org/pdf/1201.4508> の **Proposition2、printed pp.6–7** を直接読んだ。
対象はnormal・integralな曲面とample spanned H、次数H²=5であり、特異曲面を含む。
本件のH=O_S(1)はvery ampleなので同命題の全場合分けを使える。
証明中のg=1とg≥2の分岐も同頁で確認した。

S. Schröer, *Normal del Pezzo surfaces containing a nonrational singularity*,
revised 28 October 2000,
<https://www.math.uni-duesseldorf.de/~schroeer/publications_pdf/del_pezzo.pdf>
の **Theorem2.2、printed pp.4–5** を直接読んだ。
perfect field上のproper normal surfaceにnonrational singularityがあり、
全n>0でnKがeffectiveでなければq=h1(O)=0を結論する。
本件はC上であり、その射程内である。

## 2. 楕円scroll以外の全場合はq=0

Tironi Proposition2において、Hによる完全写像は本件で閉埋込みなので、
次数五のP2へのcoverの場合は除かれる。
残る楕円scroll以外の型についてq=0を示す。

有理scrollは滑らかな有理曲面なのでq=0。
有理曲線上のnormal coneは有理解消を持ち、Lerayの単射H1(O_S)→H1(O_S̃)からq=0。
P3のquintic hypersurfaceも構造完全列からH1(O)=0であり、本件の非退化性とは元々両立しない。

CI(2,3)でplaneに直接linkされた型はarithmetically Cohen–Macaulayであり、q=0。
ここでnormal surfaceはlocally CM、planeはACM、ambient CIはGorensteinなので
通常のcodimension二CI liaisonのACM保存を適用できる。
単なるcycleとしての残差ではなく、同命題が与えるscheme-theoretic linkageを使う。

Gorenstein del Pezzo型はω_S≅O_S(−H)なので、全n>0でH0(nK_S)=0。
nonrational singularityがあればSchröer Theorem2.2でq=0。
全特異点がrationalなら、normal Gorenstein rational surface singularityはDu Valであり、
最小解消f:S̃→Sはcrepantである。
従って−K_S̃=f*Hはnefかつbig、K_S̃²=5。
この滑らかなweak del Pezzo解消は有理曲面である。
例えばK_S̃がpseudo-effectiveでないのでrulednessを適用し、
正種数gのruled曲面のblowupならK²=8(1−g)−r≤0となることからg=0を得る。
Lerayにより再びq=0となる。
nonrational特異点をrationalと仮定したり、特異del Pezzoに無条件のKodaira消滅を使ったりしていない。

## 3. q=0全枝を完全Euler源で排除する

INTEGRAL_SURFACE_LINEAR_NORMALITY.mdにより、degree五のSはlinearly normalである。
q=0ならRANK_ZERO_REGULAR_LINEARNORMAL_SURFACES.mdの完全Euler源を適用できる。
degree五なので二つの独立quadricは無く、πの全liftは零又はq2 Aだけ。
実五次式は零又はq2²linearであり、既約被約Dに矛盾する。
従ってTironiの全非楕円型が排除され、Sは滑らかな楕円scrollとなる。

## 4. 全楕円scrollからscalar q52へ

Hがvery ampleかつdegree五なのでrank二bundle Fはstableである。
次数二以下のline quotientがあれば対応section上の制限がvery ampleでなくなるためである。
IRREDUCIBLE_NONNORMAL_CLEAN_RANK_ZERO_RECOGNITION.md §§7–9の
noquadric・全bundleの3regular ideal分解・Bottによるh0(Λ²T I_S)=1は
この埋込みscrollの性質だけを使い、clean/node条件は用いていない。
Markarian–Polishchukの同じscroll上で零となる非零q52がこの一次元空間を張る。
従ってπはそのscalar倍で、命題が従う。

これは既約非正規quinticの全曲面のnormality又は次数五を証明するものではない。
残りのdegree・非normal曲面・rank二restrictionは明示的な補集合として保持する。
