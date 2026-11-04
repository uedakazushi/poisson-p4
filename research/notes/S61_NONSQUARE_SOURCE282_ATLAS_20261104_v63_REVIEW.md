# source304・282 atlas と双対43：独立監査

2026-11-04。対象は frozen `S61_NONSQUARE_SOURCE282_ATLAS_20261104_v63`。
基礎は FULLU の正則整域R。新しいpivot/minor探索・GB・実239列四次行列は
計算しない。既受理の11+4関係と元六quadに結ぶ小さい疎算術だけを検査する。

## 実sourceと22列

15関係の選択source15座標は `(I11,*;0,I4)`。従ってsourceのその15座標は
関係で残る304座標へ解ける。単に対応する15行を除いた **固定304行**が、
同じimage/cokernelを与える。選択15行はgenuine cubic/quarticだけなので
最初のquadratic174行を変更しない。全index・source definitionを既受理証拠と比較する。

作者のSymPy ringやproducerを使わず、元のv-light七座標とcore三変数に直接
束縛したFraction Laurent疎辞書で、Kの全963係数を復元する。保存分母は
非零定数またはaだけ。全22関係の元六quadとの積を全係数で零と確認する。
旧cert座標とv座標の両方向linear mapの積も恒等行列である。

最初15列は元六quadのKoszul。prime六式への変換Aは単元det −80を持ち、
旧15 Koszulとの変換は ∧²A、そのdetは (−80)^5。残る六light倍と新22本目は
同じ受理済みcertをAで引き戻したものと全係数で比較する。従ってKは単なる
22個の恒等式ではなく、既受理の **whole22 kernelの同じ基底**を与える。
pointwise独立性も継承される。H-cofactorの六行は零である。

よって `E152=coker(R22→R174)` はprojective rank152。
`P282=E152⊕R130→R239` は同じcokernelを持つ。Kの22-subset minorの非零開は
FULLUを覆い、各開でその22 source行も除けば実282行表示となる。
minor自体や単一大域基底を計算済みとはしない。

大域的な自由282表示の **存在** は正しい。projective E152により左逆Pが存在し、
`R174≅E152⊕R22`。残る130列を22+108に分けた `[A,B,C]` に対し
`[A+B P,C]` がこの同型による自由282表示になる。Pは未計算である。

## targetと双対の証明

target239はpure210、H·Sym²light28、H²の直和。
quadratic174像はpure146像とH·六quadに分かれる。前者の商Q64の反射性は
独立受理済み。後者の六式の定数独立性はrank6 subbundleを与え、その商F22は
projective。従って `Q87=Q64⊕F22⊕R` は反射的、pd≤1。
残る実130行の `R130→Q87` も同じcokernelである。有限purequad例外を除き
Q87はfreeだが、全域free87とはしない。

Cを実四次cokernelとする。既受理codim Supp C≥2によりN:P282→R239は
generic全射で、双対は単射。したがって

    0→R239*→P282*→D43→0,   pd D≤1,   Ext¹(D,R)=C.

高さ1ではNは全射でDはfree。高さh≥2ではdepth D≥h−1≥1。
従ってtorsionのassociated primeはなくDはtorsion-free。
高さ3では既にdepth≥2で、高さ2でS2を満たすことはABによりfreeと同値。
ゆえに **D reflexive ⇔ Cにcurve支持がない（有限点のみ）**。
また **D locally free ⇔ C=0**。高次Extは零である。
bundle chartでFitt43(D)=I239(N)=Fitt0(C)、Fitt44(D)=I238(N)も正しい。

これらをDの反射性の証明、あるいは全射性の取得とは読まない。
Q64の反射性とD43の反射性は別である。`diag(I238,[x,y,0^42])` はcoker
R/(x,y)とabstract free43 kernelを持つ反例で、37個の零sourceを加えても
free80 kernelだけでは欠損curveを除けない。`[x,y,z,0^41]` の反例はD反射的
でも原点欠損が残る。これらは一般代数の例であって実Jacobi点ではない。

238 pivotの追加開では実cokernelは44 Schur scalarsの商となる。
rank≤237 locusはその開で覆えていない。C122 chartのSchur同型
`(p,z)↦(p+A z,C122 z)` はkernelをzへ同型に写すため、D15はpointwise独立な
rank15 subbundleとして降り、rank107 source→rank64 targetを与える。
そのcoefficient transportや固定actual107行列は未計算である。

本監査はsource atlasと条件付き代替gateを受理する。Cのcurve・点の消滅、
full239/C122 rank64、affine exhaustionは未証明のままである。
