# integral CI(2,2)の半単純pencil全体：既存symbolic証明からの拡張

2026-11-07、主義務G。S⊂P4をintegral・非退化な二quadricのcomplete intersectionとする。
pencilが非退化quadric q1を含み、q1に関するq2のself-adjoint operatorがdiagonalizableであると仮定する。
**結論。** S上で零となるrank四Poisson bivectorの実五次退化因子は、既約・被約にならない。
Sのsmoothness、normality、ordinary double又はvertical divisorの仮定は無い。

## 1. 全sourceの完全性は特異CIでも同じである

`E=Λ²T_P4` とする。q1,q2はregular sequenceなのでKoszulをEとtensorした列
`0→E(-4)→E(-2)²→E I_S→0` はexactである。
外積Euler及びline bundle cohomologyにより
`H⁰(E(-4))=H¹(E(-4))=0`、`H⁰(E(-2))=Λ²C⁵`。
従って全rank零sourceは一意に `Π=q1 A+q2 B`、A,B全constant skew matrixとなる。
これは既存smooth CI22証明のKoszul部分そのもので、Sのsmoothnessを用いない。

q1に関して異なるeigenspaceは直交し、q1のrestrictionは各eigenspaceで非退化である。
よってdiagonalizableという仮定から
`q1=Σ_i xi², q2=Σ_i λi xi²` と同時対角化出来る。

## 2. symbolic certificateが本当に反転しているもの

既存のimmutableな証明書
`reports/CLASSIFICATION_SMOOTH_CI22_ZERO_PARITY_SYMBOLIC_20261107.json`
及び別実装review
`reports/CLASSIFICATION_SMOOTH_CI22_ZERO_PARITY_REVIEW_20261107.json`
は、四index block I={0,1,2,3}の全十二parameter二次monomialについて、
`Pf(A_I), Pf(A_I,B_I), Pf(B_I)` を元の全350 Jacobi行の有理結合で表す。
保存されたrank minorは正確に

`−64(λ0−λ1)³(λ0−λ2)(λ0−λ3)(λ0−λ4)`。

保存された全functionalの分母は、非零rational定数及び
**λ0と各他のλjとの差だけ**の積である。
従って、この恒等式は他のλ同士が重複していても、λ0がsingleton eigenvalueなら有効である。
当該別実装reviewは分母を払った多項式恒等式を全五blockへ座標置換して検証済みであり、
smoothな一つの数値pencilからの推測ではない。

一般のI4と、その中に含まれるsingleton index sについて、全座標置換でsを中心index0へ送る。
同じcertificateからPoisson支持上で
`Pf(A_I)=Pf(A_I,B_I)=Pf(B_I)=0`、従って `Pf(Π_I)=0` を得る。
必要なのは中心λsの一意性であり、五つ全てのλの相異性ではない。

## 3. integral・非退化という仮定が残す全multiplicity型

pencilにrank二以下の非零quadricがあれば、それは複素数体上で二linear factorの積又は平方である。
prime I_Sがそのquadricを含むためlinear factorも含み、Sの非退化性に矛盾する。
従ってどのλのmultiplicityも三以上にはならない。
実際、multiplicity m≥3なら `q2−λq1` のrankは5−m≤2である。
q2がscalar q1となる一固有値の場合もcomplete intersectionではない。

よって半単純pencilで可能な全multiplicity partitionは

`1+1+1+1+1, 2+1+1+1, 2+2+1`

だけである。`3+1+1, 3+2, 4+1, 5` 等はここで全て排除される。

最初の二型では各I4に少なくとも一つsingleton indexがある。
従って§2により全五Pfaffianが零となり、全Poisson sourceはrank二以下である。

最後の `2+2+1` 型では、singleton indexをsとすると、
I4がsを含む四Pfaffianは零である。
残りはsを除いた一つのPfaffianだけなので、実anticanonical式は

`D_Π=(−1)^s xs Pf(Π_without_s)`

となる。これは零又はlinear factor付き五次式であり、既約・被約五次式ではない。

## 4. 閉じた射程と残す補集合

従って、特異曲面も含む全integral半単純CI22のrank零枝は、
未知G一般元の既約・被約非正規quinticから排除される。
追加の数値計算、GB又は旧証明書の書換えは用いない。
nonsemisimple self-adjoint operatorを持つpencil、及び全memberが退化するpencilはこの命題の補集合である。
`2+2+1` 型で可約退化因子を持つrank四tensorの存在・既知成分帰属を、ここからは主張しない。
