# rank二の非正規曲面：conormal Poisson moduleによる全枝排除

2026-11-07、主義務G。πをP4上の一般rank四Poisson構造、Dをその実退化因子とし、
SをDの非正規曲面の既約成分とする。S一般点でπのrankが二の場合を扱う。
ordinary node、normalization Xのsmoothness、導手の被約性又はpole次数は仮定しない。

## 1. conormalのLie bracketは零

SはPoisson不変なのでideal I=I_SはPoisson idealである。
Sの滑らかな一般点でnormal座標u,v、接座標s,tを取る。
mixed成分はIに入り、接成分π^{st}はunit。
DはS沿いに二重以上なので実PfaffianはI²に入る。
そのnormal一次項はπ^{st}とπ^{uv}のnormal一次項の積であり、
π^{uv}もI²に入る。従って `{I,I}⊂I²` が一般点で成立する。

I/I²の誘導Lie bracketはSの滑らかな部分でvector bundle間のO_S線形写像である。
一般点で零なので、S_reg全体で零となる。
この議論は横断quadratic形が非退化か退化かを区別せず、
全node・cusp・より高次の非正規横断型に適用する。

## 2. conormal bundleのflat connection

S_reg上で `N*=I/I²` とする。Sの函数fのambient liftをf̃とし、

`∇_{df}(j mod I²)={f̃,j} mod I²`

と定める。f̃の取り替えは `{I,I}⊂I²` によって影響しない。
Leibniz則によりcotangent Lie algebroidへのconnectionとなり、
Jacobi恒等式によりflatである。
S_reg上のPoissonが至る所非退化ならanchor `Ω¹→T` は同型なので、
これは通常のholomorphic flat connectionである。
従って任意の滑らかなprojective curve C⊂S_regに対して

`deg det(N|C)=0`。                                            (1)

一方、滑らかな埋込みの随伴公式はS_reg上で

`det N=ω_{S_reg}⊗O_S(5)`                                    (2)

を与える。

特にSが滑らかなcompact symplectic surfaceならω_S=O_Sなので、
(1)と(2)を一般hyperplane curveへ制限すると `0=5deg S` となって矛盾する。
従ってそのようなrank二の非正規曲面は存在しない。
これはrank零曲面の排除ではない。

## 3. normal SでresolutionのKodaira次元が非負なら全て排除

Sがnormalで、resolution S̃のKodaira次元が非負と仮定する。
normal surfaceの特異点は有限個であり、S_reg上の非零Poisson sectionは
effective Weil divisor `R~−K_S` を与える。
あるm>0でS̃のpluricanonical sectionが存在する。
resolutionはcodimension一で同型なので、そのdivisorを押し下げることで
mK_Sがeffective Weil divisorにlinear equivalentとなる。
従ってH=O_S(1)に対して

`K_S H≥0`, `R H=−K_S H≤0`。

HはampleなのでR=0であり、K_S H=0。
従ってπはS_reg全体でsymplecticとなる。
一般hyperplane curve Cは有限個の特異点を避ける滑らかなprojective curveであり、
(1),(2)から

`0=deg det(N|C)=K_S H+5H²=5deg S>0`

という矛盾を得る。

**結論。** rank二のnormal非正規曲面成分Sについて、
resolutionのKodaira次元は必ず−∞である。
SのQ-Gorenstein性、S自体のlci性又はnormalization X・conductorの仮定は不要である。
canonical Weil classと、特異点を避ける一般curveだけを使った。
非normal Sでは特異曲線を一般curveが避けられないため、この延長は主張しない。

## 4. stable楕円scrollではrank二restriction自体が無い

Sが楕円曲線B上のstable rank二degree五bundle Fによるsmooth embedded scroll
`S=P_B(F∨)`、`O_S(1)=O_P(F∨)(1)` の場合、

`H⁰(S,−K_S)=H⁰(B,Sym²F⊗det(F)⁻¹)=H⁰(B,End₀F)=0`。

最後の等号はstable bundleがsimpleであることによる。
従ってSをPoisson不変とするambient πのrestrictionは零であり、
rank二restrictionは不可能である。
この全枝排除はordinary-double又はvertical成分の条件を使わない。
既知T40の楕円scrollがrank零であることと整合する。

## 5. 射程

第3節は全normal Sの非負Kodaira次元枝の不存在、
第4節は全stable degree五楕円scrollのrank二枝の不存在という幾何学的証明である。
有限標本又は新しい計算のPASSからの外挿ではない。
normal Sのruled/rational解消、非normal S、他のrank零曲面は残る。
滑らかなrank二Sかつ一般node・vertical成分無しの全枝は、別の
`IRREDUCIBLE_NONNORMAL_RANK2_CLEAN_SURFACE.md` と全PV tangent源の独立次元証明により、
既に未知W優越成分から排除されている。
既約非正規quintic全体やG全体の完了は主張しない。
