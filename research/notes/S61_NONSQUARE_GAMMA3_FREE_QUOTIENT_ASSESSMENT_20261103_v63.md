# 全 core 上の自由な Γ₃ 商：既存結果と非循環の前提

2026-11-03。read-only の構造調査。既存証明書の再計算、GB、単位 pivot の追加、ledger の変更は行っていない。

## 結論

**現存の参照モデルから、全 `U={aFBC≠0}` 上の実際の91次元 cubic 係数空間に対する自由階数34の全射商 Γ₃ は得られなかった。** その存在が否定されたのではなく、この前提の独立な構成が未完了である。既存の generic Hilbert 関数、指定点の B65、水平自由65、自由な row module57は、いずれもこの全射商の代用にはならない。

一方、提案された行列式による条件付き論法は正しい。Γ₃ の**局所自由性**だけでも十分であり、全域で一つの基底を持つという強い意味の自由性までは必要ない。

## 実際の91列

`reports/S61_BINARY_NONSQUARE_CUBIC_REPLAY_20261103_v62/matrices.tar.gz` の member `full` が、元135式から再構成した63×91行列である。完全な hash は

```
dae59904925260605628d925ae0565de084c9af0ae1bccf1774f1af067830ad6
```

列順は同ディレクトリ `source-report.json` の top-level `normal_monomials` に保存されている（主 `report.json` では `source.normal_monomials`）。八法線変数は重み `(1,1,1,1,1,2,1,1)` を持つ。全91列は84個の light³ と7個の heavy×light であり、七変数の通常三次84列と混同していない。

係数環を `R=QQ[a,e,d,(aFBC)⁻¹]`、自由 ambient module を `V=R⁹¹` とする。`S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62` により、全63行の像 `M⊂V` は57本の明示基底を持ち、`M≅R⁵⁷` である。同定する57行の展開と六関係は全 U 上で正則であり、この結論に `dJ≠0` は必要ない。

`S61_BINARY_HYPERPLANE_CUBIC_CODIMENSION_20261102` により、繊維の階数が57未満になる支持 Σ は U 内で余次元2以上である。この定理は単なる product defect ではなく、上記の**全 cubic 行列**を対象とする。

## 条件付き行列式補題の厳密な証明

上のデータに加えて、局所自由階数34の R-module `G` と全射

\[
\pi:V\twoheadrightarrow G,\qquad \pi(M)=0
\]

が与えられたと仮定する。`N=ker π` は階数57の直和因子であり、`M≅R⁵⁷→N` は一般点で同型となる。その行列式は可逆層 `det N` の非零切断である。

高さ1の点では Σ が存在せず、`M→V` の繊維階数は57である。従って `M→N` の行列式は全ての高さ1点で単位となる。正則整域上の可逆層の非零切断がどこかで消えるなら、その点で局所的に表す非零非単位元は高さ1の素因子を持つ（主 ideal 定理）。よって切断はどこでも消えず、`M→N` は全 U 上同型となる。従って `V/M≅G`、全点で cubic rank57となる。

`G` が自由なら `N` は projective な直和因子であって、自由と先に仮定する必要はない。上の局所的な行列式の議論がそのまま適用できる。

逆に全点で rank57が既知なら、非零57 minor の開被覆で通常の pivot を行うことにより `V/M` は局所自由34となる。従って、既存の余次元2定理の下で、このような局所自由な全射商を**存在だけ仮定すること**は、実質的に未解決の全点 rank57を別表現で仮定することになる。独立な幾何構成や有限の全射証明書が必要である。

## 既存候補の範囲

|既存資料|実際に証明された内容|全 U の Γ₃ に足りない点|
|---|---|---|
|`S61_BINARY_TRANSVERSE_ALGEBRA_20261101`|指定された固定β点の横断 Artin 商 B65、通常 HF `(1,7,23,34)`|一点の横断商であり、全 core に渡る商写像・平坦性はない|
|`S61_BINARY_FORMAL_RING_20261101` と review §4–5|指定点と正当な一定階数開で完成環が3変数底上自由65、`I³` が自由34|特殊横断長65の上界を使って初めて自由性を証明する。全 U にこの上界を与えていない|
|`S61_BINARY_UNIFORM_CUBIC_20261101` §4–5|同じ一定階数 incidence 上の通常三次関係の上界50、法線長の下界65|関係の上界は quotient の生成元34本・全射・平坦性を与えない。84通常三次と91重み三次も直接同一でない|
|`S61_BINARY_HYPERPLANE_WEIGHTED_CHART_20261102`|全 U 上の正確な八変数 weighted implicit model|モデルの正則な座標表示であって、低次商の自由性を結論していない|
|`S61_BINARY_HYPERPLANE_GENERIC_RANKS_20261102`|一般 weighted HF `(1,7,23,34)`、全点の cubic 上界57|全点下界57や全点自由34は明示的に未主張|
|`S61_BINARY_HYPERPLANE_BFREE_HORIZONTAL_20261102` §5 と review|`K=QQ(a,e)` 上、全 d の weighted 代数が `K[d]` 自由65、次数3は自由34|`a,e` の関数体化で縦曲線の支持を除いている。全 U の降下ではない。低次階数に既存 codimension theorem を用いる|
|`S61_BINARY_DEEPER_CORE_OPEN_FORMAL_20261101`|別の深部七変数 source の指定 core 開で相対自由65|`F=0` に位置する深部 source の定理であり、現在の `F≠0` 三 core 全域のモデルではない|
|`S61_BINARY_NONSQUARE_CUBIC_ROW_BASIS_20261103_v62`|全 U の**row module**が自由57|埋込みの特殊化で階数が落ちないこと、余核が局所自由であることは含まない|

滑らかな実在 binary incidence は被約な台と三つの tangent 切断を与える。しかし、その滑らかさ自体は非被約な法線厚化の次数3を自由34にしない。ここで求める Γ₃ と、その滑らかな被約 incidence を区別する。

## 同じ数値を持つ最小反例

一般階数、自由 row module、余次元2の defect だけでは不足する。`A=QQ[x,y,z]` 上の写像

\[
A\longrightarrow A^2,\quad 1\longmapsto(-y,x)
\]

は単射で、その像は自由1である。一般階数1、defect は余次元2の直線 `x=y=0` である。その余核は ideal `(x,y)` と同型であり、自由1ではない。自由な候補 target `A` に向かう写像 `(u,v)↦xu+yv` は上の関係を消すが、像は `(x,y)` に過ぎず、全射ではない。

56本の identity block と33本の自由な余核を直和すれば、**自由57→自由91、一般余核34、余次元2 defect**を持つ同じ大きさの例になる。この例の余核は torsion-free であり、double dual は自由34だが、元91から double dual への写像は defect 上で全射でない。従って saturation、torsion-free 化、reflexive hull の語だけで欠けた全射性を補うこともできない。この例は現在の Poisson 行列の欠損を示すものではなく、推論の前提を区別するための反例である。

## 独立構成に必要な有限証明書

最も直接的な新しい候補は、実際の91列順で `P:R⁹¹→R³⁴` と右逆 `T:R³⁴→R⁹¹` を与え、

\[
P M=0,\qquad PT=I_{34}
\]

を QQ の全係数で確認するものである。全分母は既知 `a,F,B,C` の積だけでなければならない。局所 chart の有限被覆で同じことを行い、target が局所自由34であることを別に証明してもよい。一般点の kernel basis の分母を払うだけでは `PT=I` を失い得る。

幾何的な代替は、実際の weighted normal algebra を商とする別の厚化 Γ を構成し、次数3の34個の基底と正則な reduction rule を、現在の full cubic rank57やその長65上界に依存せずに証明することである。その際、三つの実在 binary tangent 切断との整合、heavy 座標の重み2、135元式に由来する全63行の消滅、core 全域の平坦性が必要となる。

この調査で、そのような `P,T` 又は独立な全 U の厚化模型は発見できなかった。現在の実行済み単位還元を超える rank 下界、J=0 defect の有限性、quartic 被覆、ambient 飽和除去は結論しない。

## 保存した read-only receipt

`code/s61_nonsquare_gamma3_assessment_20261103_v63.py` は全91単項式の重みと完備性、保存 full member の hash、既存 report の適用範囲、上記最小反例を軽量に検査し、調査した依存の hash を記録する。既存の135式再構成や GB は再実行しない。保存先は `reports/S61_NONSQUARE_GAMMA3_FREE_QUOTIENT_ASSESSMENT_20261103_v63.json`。
