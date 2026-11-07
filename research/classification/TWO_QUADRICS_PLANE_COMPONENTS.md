# 二つの平面conic配置から得る二つの23次元全Poisson成分

2026-11-07、課題E_reduced_nonnormal。
`TWO_QUADRICS_PLANE_CONIC_CONFIGURATIONS.md` の完全な七配置・全六磁気源を入力とする。
本稿は以下の二源のPGL飽和像閉包が、全Poisson空間の相異なる23次元成分であることを示す。
別の三源の横断支持は非孤立であり、その一源は後続の全β=0座標変更で被覆した。
残二源の成分帰属は未解決。精密なnilpotent厚さや一般軌道次元は主張しない。

元の平面座標をX,Y,Z、`q=XZ−Y²` とする。

| ID | H | q₁,q₂ | 幾何 |
|---|---|---|---|
| contact31_common_bitangent | X+8Y+16Z | q,q+XY | conic交点3+1、Hは両conicに異なる点で接する |
| contact4_other_tangent | Z | q,q+X² | conic交点4が一点、Hはその点外で一方に接する |

両源は全六係数の磁気二次式を持つ実Poisson族であり、全係数のstrict Jacobi・実Pfaffianは
既に別Fraction実装で確認済みである。各源の全imageが一つの成分へ含まれることも
以下のsource閉包認識から得る。単一代表の軌道閉包だけの成分性ではない。

## 1. 全150係数での実次元下界と真正slice

五変数の全二次bivector係数150個を取り、全projective Jacobi350係数を用いる。
h-Casimir又はstrict Jacobiの条件を周囲の未知変数へ追加しない。
全源の六磁気係数を(1,2,3,5,7,11)へ置いた点をwitnessとする。
非零係数の有理分母を共通整数scaleで払っただけであり、源を小さくしていない。

全Euler-vertical24次元、全GL5作用、全体scaleと全六の実parameter方向を保持する。
これらの行列の正確なrankは、両witnessで

`vertical=24`, `vertical+group+scale=47`, `full_actual_source=48`。

したがって全sourceのaffine lift像閉包は次元少なくとも48、
vertical24とscale1を除いた射影像閉包は次元少なくとも23である。
全方向は実際の多項式族又は群作用の微分であり、単なる形式一次方向ではない。

Jacobi微分をD、実source接空間をFと書く。
両例ともrank D=101、rank F=48である。
kernelを全て計算し、literal normal Rを選ぶと

`ker D=F⊕C R`、`dim ker D=49`。

Dの独立な101列が指定する標準座標方向とRは102次元の実線形部分空間Sを張る。
Dは前者に単射、FとRを零に写すので `C¹⁵⁰=F⊕S` である。
Sのannihilatorを取った48個の線形cutが真正なsource横断sliceを与える。
選択した101×101の非零minorにより、その101座標変数を形式的に消去できる。
残る一変数をtとする。

## 2. common-bitangent源の二次支持

`CLASSIFICATION_TWO_QUADRICS_PLANE_NORMALS_20261107.json` に
全150列を零にするliteral左cokernel汎関数ℓを保存した。
選んだRについて

`ℓ D=0`、`ℓ Jac_projective(R)=−120295458720≠0`。

従って真正sliceの残余方程式は非零係数のt²を初項に持つ。
その被約支持は原点であり、sliceの次元は零である。

独立 `--plane-components` checkerはproducer/CASをimportせず、
元150係数から全D列を別に作り、全六実parameter方向、全vertical/group方向、
完全kernelの張り、101minor、ℓの全150列への零作用と上のpairingを確認した。
報告は `CLASSIFICATION_TWO_QUADRICS_PLANE_COMPONENT_REVIEW_20261107.json`。

## 3. other-tangent源の三次支持

この源の一法線二次cokernel障害は零である。
そこで、同じ真正sliceで、選択した101 Jacobi行を零にする**固定された**二次補正Tを
101個のliteral標準座標方向内で一意に求めた。
Tは既知source方向や任意のkernel補正を加えたものではない。
`P+tR+t²T` の三次Jacobi係数をCとする。
実際のcokernel汎関数ℓについて

`ℓD=0`、`ℓ Jac_projective(R)=0`、`ℓC≠0`

が成り立つ。これにより、sliceの一つの残余方程式は非零t³を初項に持つ。
選択行を三次でさらに補正しても、その項はℓD=0によりpairingを変えない。
従って真正sliceの被約支持はここでも原点である。
全kernelに対する任意補正の不変性や形式的積分解は、この証明の追加仮定ではない。

`CLASSIFICATION_TWO_QUADRICS_PLANE_CUBIC_SUPPORT_20261107.json` は
固定T、全三次係数ベクトル、ℓの全係数とpairingを保持する。
別Fraction実装は、元Dの150列・選択minor・Tの指定座標内性・二次補正恒等式・
全三次係数と全cokernel作用を再構成し、非零の実t³係数を確認した。
結果は `CLASSIFICATION_TWO_QUADRICS_PLANE_CUBIC_REVIEW_20261107.json`。
任意の四次以上のjetやexact nil環の同定は行っていない。

## 4. 全空間成分への帰結

各真正48-cut sliceは零次元なので、Krullの不等式から
元の全150係数の周囲局所次元は48以下である。
実source像の下界48と一致する。
実sourceのPGL・vertical・scale飽和は既約であり、その像閉包はwitnessを含むので
次元48の全affine成分そのものである。
vertical24とscale1を除くと、射影Poisson空間の**23次元全既約成分**となる。
よって各完全六磁気源の全点がその成分に入る。

ここで一般非被約環のexact同型、一般軌道の稠密性、witnessでの局所成分の一意性は
結論していない。二次又は三次の横断初項だけで必要な被約成分認識が得られる。

## 5. 既知成分・相互との相違

両sourceの全退化五次式は本質変数数≤3を持つ。
五偏微分の係数行列rank≤3は閉条件であり、PGL飽和閉包でも保たれる。
非零の五次式が三変数だけに依存すると、P4内のvertex P¹で重複度5を持つ。
四次元のSNC divisorの重複度は高々4なので、この閉条件内にはSNC点がない。
従ってSNC点を持つMPS40のどの成分とも一致しない。

一般rank四により六rank二成分と異なる。
射影次元23により、次元22,24,25,28の既知追加成分と異なる。
とくにCPLUS_EXTは次元24なので、一般同因子型でも次元だけで区別できる。
同次元23の既知非MPS成分の一般退化型は、既存成分監査により
1+4、五超平面又は非被約、1+1+1+2、2+3、1+1+3であり、
新しい1+2+2型とは異なる。
入口は `KNOWN_COMPONENT_AUDIT.md`、`ADDITIONAL_RANK4_LATE_COMPONENT_AUDIT.md` と
そこから参照される現行個別台帳である。

新しい二つ同士は、一般Dを共通vertexで商した平面conicの交差数分割が3+1と4で異なる。
因子を交換してもこの分割は変わらない。
各sourceのdense constructible像にはこの配置型を持つ非空開があるため、
両閉包が同じ既約成分ならその二開が交わることになり矛盾する。

## 6. 三つの非孤立結果と、その後の被覆

以下は同じ全150係数・全六磁気方向を保ったまま、二次と一回の必要な三次支持を計算した。
独立Fraction実装でも同じ非孤立結果を確認した。

| 源 | 二次初期式 | 三次の後の残支持 |
|---|---|---|
| contact31_tangent_at_simple_intersection | t₀²−625t₁², t₀t₁−25t₁² | t₀=25t₁が残る |
| contact4_secant_at_contact | なし、一法線 | 三次も零で法線直線が残る |
| contact31_chord | なし、二法線 | 三次二式の正確なgcdはt₀+t₁、直線t₀=−t₁が残る |

これらのliteral法線軸、及び第一行の唯一の二次支持直線を調べても、
Jacobiを満たす直線Poisson pencilにはなっていない。
従って全実変形の存在も、成分認識もまだ結論しない。
このjet経路は三次で停止した。
後続 `TWO_QUADRICS_PLANE_EVEN_SUPPORT.md` はcontact4_secant_at_contactを
全六磁気係数でβ=0へ送る正則座標変更を構成し、既知有限閉和への全source所属を証明した。
残二源は全sourceの構造的変形又は既知成分認識を行う。
非平面1+2+2と他因子型も未解決のまま保持する。
