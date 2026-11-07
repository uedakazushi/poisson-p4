# normal degree五rank零曲面の認識：別作者による適用範囲監査

2026-11-07。producerはroot、reviewerはgeneric_exhaustion。
対象は `NORMAL_DEGREE_FIVE_RANK_ZERO_SURFACE.md` の幾何証明である。
**判定：記載された射程でPASS。** normal・integral・非退化なdegree五rank零曲面Sがあれば、
既約被約実quinticのrank四Poisson構造はT40となる。
Sのsmoothness又は横断node性、正規化/導手のsmoothness、vertical divisor不存在を要しない。
以下はexact CASによる再現ではなく、一次資料と仮定を照合した幾何学的な別読解である。

## 1. normal曲面全体を扱う分類入力

A. L. Tironi, *Normal projective varieties of degree 5*,
arXiv:1201.4508v1, <https://arxiv.org/pdf/1201.4508> を直接読んだ。
Proposition2（printed pp.6–7、PDF index5–6）は、任意の代数閉体上の
irreducible reduced normal surfaceとample spanned H、H²=5を仮定する。
本件のHはvery ampleなので全仮定を満たす。
statementの七場合と、証明のsectional genus0/1/≥2の分岐を実際に照合した。

このsurface propositionを使えば、導入部の高次元定理に現れる
generalized coneの意味を別途補って分類を拡張する必要が無い。
例えばgenus二degree五曲線のconeが現れても、normalな場合は同命題のCI(2,3)link型に含まれ、
「coneなら有理解消」と誤って扱っていない。
有理解消を使うのは、statementの**smooth rational curve上のcone**という明記された別場合だけである。

Hによる完全写像はこのSの閉埋込みなので、degree五P2 coverは除かれる。
有理scroll、smooth rational curve上のcone及びquintic hypersurfaceはq=0である。
CI(2,3)でplaneに直接linkされた型も、normalゆえlocally CMで、
codimension二CI liaisonのACM保存によりq=0となる。
これはdegreeだけから残余cycleをlinkと同一視した議論ではなく、同命題が与えるscheme-theoretic linkageを使う。

## 2. Gorenstein del Pezzoの場合のq=0

Tironiの場合5は `ω_S≅O_S(-H)` と明記しており、SはGorensteinである。
従って全n>0でH⁰(nK_S)=0、またnormal surfaceの双対性からp_g(S)=H²(O_S)=0。

nonrational singularityが一つでもある場合、S. Schröer,
*Normal del Pezzo surfaces containing a nonrational singularity*, revised28 October2000,
<https://www.math.uni-duesseldorf.de/~schroeer/publications_pdf/del_pezzo.pdf>
のTheorem2.2（printed pp.4–5、PDF index3–4）を独立に直接読んだ。
同定理はproper normal algebraic surface、perfect ground field、
nonrational singularity及び全n>0でnKがineffectiveという仮定だけで、q=0を結論する。
定理にminimal resolutionを導入する記述はあるが、SのPicard rank一又は全Sのsmoothnessを仮定しない。
著者のqはPic0の次元であり、p_g=0の下でh¹(O)に一致することも同論文§1で確認した。
従って本件への適用は有効である。

全singularityがrationalの場合は、Gorenstein rational surface singularityがDu Valであるという
標準的なsurface singularity定理を使う。
minimal resolutionはcrepantで、`−K_S̃=f*H` はnefかつbig、K²=5。
Kがpseudo-effectiveでないためS̃はruled又はrationalであり、
正種数gのruled surfaceのblowupなら `K²=8(1−g)−r≤0` に反する。
従ってS̃はrational、q(S)=0である。

Tironi Proposition2のstatement単独からdel Pezzoのq=0又はh⁰(H)=6を仮定してはならない。
現証明では上記Schröer入力を明記することでこの橋を確保している。
同入力を取り除くには、normal Gorenstein del Pezzoの別のvanishing又は分類証明が必要であり、
smooth Kodaira vanishingを特異Sへそのまま適用して省略していない。

## 3. q=0の排除と、楕円scrollの全埋込み

degree五はprojected Veronese例外ではないので、別途一次資料確認済みの
`INTEGRAL_SURFACE_LINEAR_NORMALITY.md` によりlinearly normalである。
q=0の場合、Euler完全源は零又は `q₂A` のみとなり、実Dは零又は `q₂²linear`。
従って既約被約quinticに反し、非楕円型は全て排除される。

残ったTironiの場合2はsmooth P1 bundle over elliptic curveである。
Hのfiber次数一からrank二bundle Fのtautological embeddingとなり、degF=5。
Fが非stableなら次数二以下のline quotientを持ち、そのsectionへのHの制限がvery ampleでなくなる。
従ってFはstableである。

既受理の `IRREDUCIBLE_NONNORMAL_CLEAN_RANK_ZERO_RECOGNITION.md` §§7–9を再読した。
そのnoquadric、3-regularity及び全bundleのresolution
`0→O(-5)→O(-4)⁵→O(-3)⁵→I_S→0` はscrollとembeddingだけの主張で、clean/node仮定を使わない。
`E=Ω²(5)` とtensorするとBottによりΩ²(1)、Ω²(2)はacyclic、H²(Ω²)=Cなので
`h⁰(E I_S)≤1` を得る。
Markarian–Polishchukの任意のstable rank二degree五bundleに対する同じscroll上で零となる
非零q52が存在するため、全section空間はちょうどその一次元である。
従って今回のπはscalar q52であり、T40認識にclean pole議論を持ち込む必要は無い。

## 4. 保持する補集合

この定理はnormalなdegree五rank零曲面を全て扱うが、
degree六、非normal曲面、rank二restriction及び他のdegreeの曲面を分類しない。
Sが複数ある場合も、その一つが仮定を満たせば同じ認識が成立するという部分定理である。
既約非正規quintic全体又はG全体の完了を意味しない。
