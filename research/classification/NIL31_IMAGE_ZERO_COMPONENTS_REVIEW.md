# nil31 image内二源：相異なる二つの25次元成分の独立確認

2026-11-07。課題 `E_multiple_hyperplane`。
`NIL31_IMAGE_ZERO_COMPONENTS.md` と `NIL31_IMAGE_THIRD_COMPONENT.md` の成分性・相違を受理する。
half は滑らかな25次元点で証明し、third は真正二法線の零次元支持で証明した。
third の滑らかさ・局所成分一意性・精密な非被約構造は主張しない。
既存92成分に対して相異なる二成分を与えるという相対的結論であり、全分類の完了ではない。

## 1. 全源と旧模型への橋

両源について、六個の自由磁気係数 A0,…,A5 と u,v,t、最後の p 又は q の全十変数を保持した。
旧15変数源の各係数を独立に homogenize し、パラメータを

    (n0,n1,n2,n3,n4,n8,n12,n13,n14)
       = (A3,A0,A2,A5,A1,A4,u,v,t)

へ移した。half は n10=−1/2,n5=p,n6=n7=n9=n11=0、
third は n10=−1/3,n6=q,n5=n7=n9=n11=0 である。
射影置換 `(h,a,b,c,d)=(old_z,old_h,old_x,old_y,old_w)` の行列式−1と、
Euler gauge `E∧(0,0,−rc,−rh,0)` の減算を全150係数で再構成した。
非線形座標変更を PGL 同値として使っていない。

独立に得る全源は

    A=A0h²+A1ah+A2dh+A3a²+A4ad+A5d²,
    K=3r+t+1, S=ua+vh+(t−r)d,
    Pab=rac, Pac=rha, Pad=0,
    Pbc=A+K(hb−c²/2), Pbd=g+cS, Pcd=hS,

half では g=pa²、third では g=qha。
全十 strict Jacobi が全パラメータで恒等零であり、

    Pf_half=(p/2)ha³, D_half=(p/2)h²a³,
    Pf_third=(q/3)h²a², D_third=(q/3)h³a²

を文字通り再展開した。
従って p≠0、q≠0 はそれぞれ一般階数四の条件であり、他の係数条件は付けていない。
零境界も同じ既約な十変数多項式源に残る。

固定点 `A=(1,2,3,5,7,11),u=2,v=3,t=5,p/q=13`、tensor12倍において、
全十個のパラメータ偏微分を**源の式から**別に生成して producer の全係数と照合した。
補助の κ や r の微分を actual F の下界へ加えていない。

## 2. 独立算術と全周囲行列

実装は `code/classification_nil31_image_zero_components_review.py`、third は `--third`。
producer、producer helper、CAS、NumPy、legacy arithmetic は import しない。
全源は独立な LaurentRing の Fraction 疎演算で生成し、
全350式と全150方向は別作者の nil22 independent reviewer に由来する疎係数演算から再生成した。

D は全350 projective Jacobi 式の全150係数微分、V は25列の Euler 鉛直方向、
G は V と全25 GL5 方向及び scalar、F は G と全十実パラメータ方向である。

| 検証 | half | third |
|---|---:|---:|
| rank D | 100 | 98 |
| rank V | 24 | 24 |
| rank G | 47 | 48 |
| rank F | 50 | 50 |
| ker D 次元 | 50 | 52 |
| 射影接次元 | 25 | 27 |

全四行列の全係数 hash と、保存された最大 QQ minor の正確な行列式が一致した。
全61列の actual F に対する `DF=0` も直接確認した。
half では非零100 minor と50個の独立な kernel 方向が上下から rank D=100 を確定する。
source saturation の射影次元は少なくとも50−24−1=25、全接空間はちょうど25。
従ってこの点の局所環は正則で、唯一の25次元成分を通る。
既約な全源像閉包はその成分と同次元なので等しい。
全源が singular な CA/CB 特殊点に接続するという事実から包含を推論していない。

## 3. third の真の線形 cut と二次初期支持

third の二法線は全 `ker D / im F` に補であり、全 `DN=0`、
`rank[F N0 N1]=52` を独立に確認した。
D の98個の独立標準座標方向と二法線を合わせた L は100次元。
全150係数空間で `rank[F L]=150` を確認したので、L は実際の50本の独立線形式の共通零点である。
固定因子源の一切片に任意の二法線を置いただけの判定ではない。

保存された98次の implicit D minor の正確な非零行列式を再計算した。
二本の正規化左汎関数について、**全150列**で `ℓD=0`、
法線二単項式の全三 pairing が次の係数列に一致することを確認した。

    [1,0,−113435/3825851],
    [0,1,214061/7651702].

従って真の cut での残余式の二次初期式は

    q1=s²−A t², q2=st+B t²,
    A=113435/3825851, B=214061/7651702.

Fraction 演算で

    B²−A=−1690119521019/58548543496804≠0,
    t*q1−(s−Bt)*q2=(B²−A)t³

を検証した。その根基は(s,t)である。
形式陰関数の二次補正は D の像に入るため、全 D を殺す二汎関数の初期式を消せない。
よって真の100次元 cut の Poisson 支持は零次元。
Krull の次元不等式により全150係数模型の局所次元は50以下となり、actual F50 の下界と一致する。
既約な全源 PGL 飽和像閉包は25次元の全空間成分である。

ここでは source 方向を固定した真の cut の次元だけを用いる。
source×normal の形式的積構造や、third の局所成分一意性は結論しない。
接次元は27のままであり、滑らかさを主張していない。

## 4. 既存92成分との相違

監査時点の全92成分の ID、次元、一般階数、SNC 条件、一般因子と証拠参照を
`CLASSIFICATION_NIL31_IMAGE_ZERO_PRIOR92_20261107.json` に固定した。
全92行の次元に欠落はなく、同次元25の成分は正確に

    T01,T39,T40,positive-weight NEW25,CN31独立(2,2,1)

の五つ。他の87成分は、両新源の全空間成分性を確定した後には次元で区別できる。
一般 rank 二成分とは閉 rank 条件でも区別する。

両新源の rank 四部分は、全パラメータ・全 PGL 飽和で `D=H³L²` を満たす。
この条件は proper な射影積写像

    P(V*)×P(V*) -> P(Sym⁵V*), (H,L) -> H³L²

の像であり閉じる。零の退化式も含めた affine cone は閉なので、
源の**成分閉包全体**にも保たれる。
像の非零元は全て非被約であり、SNC の reduced divisor を持たない。
従って SNC locus と交わる既知40成分、特に同次元の T01/T39/T40 とは一致しない。
成分性を確定した後の全閉包比較であり、ある特殊点が MPS に入らないという推論ではない。

positive-weight NEW25 の一般因子は3Hと滑らかな rank5 quadric。
一意分解性から、既約 quadric を平方の一次式へ変える `H³L²` 表示は不可能である。
CN31 の一般因子は三つの独立超平面で重複度2,2,1、三重一次因子を持たない。
従って五つの同次元成分も全て区別でき、両新成分は既存92のいずれとも相異なる。

## 5. half と third の相互識別

一般因子の次数・重複度は両者とも3,2なので、それだけでは区別できない。
唯一の三重超平面を使う。
与えられた projective Poisson class と Poisson 超平面 H に対して、
H-Casimir 二次代表は一意である。
実際、二代表の差を `E∧Bx` とすると H=x0 の座標では
`x0(Bx)_j−x_j(Bx)_0=0`。これより B は scalar で、差は零となる。
従ってその divergence field は H と π から intrinsic に復元され、PGL と非零 scalar に自然である。

全パラメータで、置換後の共通の projective cofactor 方向は

    N=c∂b+h∂c

となることを、gauge 前・後の divergence から独立に確認した。
half の唯一の三重超平面は a=0（old_h=0）で、その Casimir 代表の modular 場は N。
third の三重超平面は h=0（old_z=0）で、その Casimir 代表の場は `−(2/3)N`。
行列として `rank N=2,rank N²=1,N³=0` である。

half では covector a が `im Nᵀ` に入らず、a=1 上に N の零点を持つ。
third では `h∈im(Nᵀ)²` で、N の零点は h=0 に全て含まれ、h=1 で零点を持たない。
この全源の dual flag は Fraction の rank 比較でも再構成した。
したがって唯一の三重超平面を保つ affine/PGL 同値で両源は交わらない。

もし両25次元成分が同じなら、各非零 p,q の PGL 飽和 constructible 密像は
その成分の非空開を含み、二つの開は交わる。
そこで intrinsic な三重超平面と cofactor flag が上の二条件を同時に満たすことになり矛盾する。
よって half と third も相異なる。
固定 rank 上の image flag が全退化に対して閉じるとは仮定していない。

## 6. 証明書と検証状態

half の独立証明書は
`CLASSIFICATION_NIL31_IMAGE_ZERO_COMPONENTS_REVIEW_20261107.json`、
契約は `jobs/nil31_image_zero_components_review.json`。実測0.433秒、最大RSS 18,432 KiB。
third は
`CLASSIFICATION_NIL31_IMAGE_THIRD_COMPONENT_REVIEW_20261107.json`、
契約は `jobs/nil31_image_third_component_review.json`。実測1.766秒、最大RSS 21,136 KiB。

旧半源 reviewer の bytes は元の独立 report に完全保存し、共通 checker へ `--third` を追加した。
immutable producer・既存証明書は変更していない。
上記成分性・全支持包含・相互相違を受理する。
外部文献上の新規性、全成分数の上界、全分類、精密な非被約構造は未確認・未完成である。
