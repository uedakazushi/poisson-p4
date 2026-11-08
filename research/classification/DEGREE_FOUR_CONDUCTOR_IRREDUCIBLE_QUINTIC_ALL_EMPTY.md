# integral degree四導手曲面はactual既約quinticに存在しない

2026-11-08、PRIMARY E。全支持分類と受理済みpointwise排除の純粋aggregate。
PVの別作者module／Borel receiptも含め、未受理の計算枝は入力に残さない。
新しい計算、GB、minor、jet又はsource probeは無い。

**定理。** P4の一般rank四Poisson構造の実退化因子Dが被約既約quinticなら、
Dの非正規曲面支持にintegral degree四成分Sは無い。
Sのsmoothness／normality、restriction rank零／二、normal node、Rの被約性、
clean、vertical無し又はnormalization-flatを仮定しない。

## 1. 共通入口と退化したspan

Hamiltonian導分が正規化と導手を保つため、SはPoisson不変。
その一般restriction rankは零又は二である。
S⊂SingDなのでactual五次FはSのsymbolic squareに入る。

Sがhyperplaneに入るなら、prime CI idealはI=(ell,f4)。
I²はI-primaryなのでI^(2)=I²、degree五よりF=ell²A3+ell f4*constant。
従ってell因子があり、Dの仮定に反する。この議論は両restriction rankに有効。
以下はSがP4をspanする場合だけを扱う。

## 2. 全integral非退化degree四surfaceの分類を戻す

`DEGREE_FOUR_RANK_ZERO_SURFACE_EXHAUSTION.md` §1–5のclassification部分はrank非依存である。
外部almost-minimal-degree theoremの全integral scopeにより、SはACM del Pezzo、
又はminimal-degree surface Y⊂P5の外点からの単純射影である。
normality／linearnormality／generic projection centreを追加仮定しない。
ACMの場合はschemeとしてCI22。
残りのYはVeronese、S13、S22、S04の全四型である。

| 全centre／surface branch | 帰着又はpointwise結論 |
|---|---|
| normal・linearlynormal ACM、各projectionのACM centre | CI22 |
| Veronese、centreはclosed Sec内の外点 | CI22 |
| Veronese、closed Sec外 | 唯一のsmooth projected Veronese orbit |
| S13/S22、非ACM全四centre orbit | rank四quadric上の非CI quartic、actual quadric因子 |
| S04、closed Sec(C4)内のcentre | ACMなのでCI22 |
| S04、closed Sec(C4)外の全連続centre | smooth rational quartic cone、rank四quadric上の非CI surface |

surface自身のcentreは射影定義の外であり除く。
まだintegral・degree四・非退化を保つ全collision／special centreも表のどこかへ戻る。
generic centreのclosureだけで境界を覆ったとの推論ではない。

## 3. rank二の全三つの残余

**CI22。** `INTEGRAL_CI22_RANKTWO_IRREDUCIBLE_QUINTIC_EXCLUSION.md` の全normal／nonnormal theoremを使う。
nonnormalのregular22_1、nil32、singularM1、cone-node、cone-cuspの全五source、
全Jacobi／actualD／zero補集合を戻し、全actual既約quintic枝を空としている。

**非CI scroll／cone。**
`QUARTIC_SURFACES_ON_QUADRIC_DOUBLE_QUINTIC_EXCLUSION.md` と
`QUARTIC_QUADRIC_DOUBLE_QUINTIC_SEPARATE_REVIEW.md` を使う。
表の全非CI配置はrank四quadric Qを持つ。SのclassはCl(Q)=Z²で(1,3)又は(3,1)。
QがFを割らないなら、D|Q−2Sは(3,−1)又は(−1,3)となりeffectiveでない。
ゆえにQがFを割り、実既約Dに矛盾する。
S04の全outside-secant連続族でもcontaining quadricのrank四を確認済みであり、
rank三quadricへこのWeil-class議論を流用しない。

**smooth projected Veronese。**
`PROJECTED_VERONESE_TANGENT_SL2_EXHAUSTION.md` と
`PROJECTED_VERONESE_TANGENT_SL2_SEPARATE_REVIEW.md` の全十五源を使う。
全Poisson cone次元二以下と実SL2作用から、全射影点はBorel eigenlineへ入る。
完全Sym2／Sym4／Sym6分解の全三branchはrank二以下、既済rank零、actual超平面因子となり、
全rank四・actual被約既約quintic枝は空。
新しい点を選んだ検査やdimension≤22だけをEへ読み替えたものではない。

これでrank二の全degree四branchを尽くす。

## 4. rank零と全Eへのscope

rank零については既受理 `E_DEGREE_FOUR_RANKZERO_IRREDUCIBLE_QUINTIC_ALL_EMPTY` と
`CI22_NIL32_AND_DEGREE_FOUR_RANKZERO_ROOT_REVIEW.md` を別入力とする。
rank零nil32六源をrank二十二源へ読み替えず、両方の独立した全source排除を用いる。
第1節の二つのrestriction rankがともに排除され、定理が従う。

結論は特殊βを含む**全点**である。既存のG一般β排除からの外挿ではない。
degree一・二・三・五・六の導手支持、他のactual factor型、全E又は総component数は未主張。
新componentの増分は零。
全依存のexact SHAは同名manifestに固定し、旧source、証明書、主稿、global台帳、Gitは変更しない。
