# H＋二既約quadric：純冪零cofactor十型からのG用還元

2026-11-07。課題G。λ=0・joint半単純rank零ではZ=X₁+X₂もnilpotentである。
最高quartic Pf(Q)≠0から線形部A=0は不可能なので、Zのaffine共役型は
`HYPERPLANE_QUARTIC_NILPOTENT_SOURCES_REVIEW.md` の全十型である。
その全90係数源・全Jacobi・実Pfの独立再構成を共通入力とする。
ここで対象は相異なる二既約quadric q₁,q₂の積であり、quartic既約性を仮定しない。

## 1. 全十型の射程

| 型 | cofactor | H+q₁+q₂での結論 | 根拠 |
|---|---|---|---|
| 1 | Jordan4,c=0 | 全PGL像次元≤24 | 第2節 |
| 2 | Jordan4,top4 | 空 | 全Jacobi支持でPfは零又は一quadricの平方 |
| 3 | Jordan31,c=0 | 空 | Jacobi前の全Pfがzで割れる |
| 4 | Jordan31,top1 | 全像次元≤23 | `HYPERPLANE_QUARTIC_NIL31_SHORT_SUPPORT.md` と独立V2 review |
| 5 | Jordan31,top3 | 空 | 第3節 |
| 6 | Jordan22,c=0 | 全て平面七配置、像次元≤24 | `TWO_QUADRICS_NIL22_BASE_SUPPORT.md` |
| 7 | Jordan22,top2 | 全像次元≤24 | 全20源・二段fiber・全criticalcoreの独立review完了 |
| 8 | Jordan211,c=0 | 空 | 第4節 |
| 9 | Jordan211,top1 | 全て平面七配置、像次元≤24 | 第5節 |
| 10 | Jordan211,top2 | 全PGL像次元≤24 | 第6節 |

型2の結論は `HYPERPLANE_QUARTIC_NIL_REGULAR_SUPPORT.md` 第1節の
全二Jacobi分岐をそのまま使う。型3のz因子は独立全source reviewerが無条件に確認した。
型4の次元定理の全Jacobi支持と例外minor処理は別証明であり、本稿で再計算したとはしない。
型7は `HYPERPLANE_QUARTIC_NIL22_LONG_SUPPORT_REVIEW.md` により、全Jacobi・
全670lift・全critical層・一様fiber・affine軌道を別作者が独立再構成して閉じた。
従って全pure-nil十型が未知W優越一般成分から排除された。
全E成分帰属又は他因子型のG完了をこの表から主張しない。

## 2. 型1：残る二quadric枝と弱い次元上界

`HYPERPLANE_QUARTIC_NIL_REGULAR_SUPPORT.md` 第2–4節の全scalar分岐、
実際のN²-centralizer変換、全閉補集合の処理を入力とする。
その還元ではquartic既約性を使う前に全実Pfaffianを列挙している。
最高quartic零、線形因子、quadric平方、絶対既約quarticの各枝を除けば、
相異なる二既約quadricが残り得るのは

`n₁₂=−1/2`, `n₁₀=0`,
`Pf P=I₂²/48−n₀w²/2`, `I₂=2wy−z²`。

ここではn₀≠0である（零なら同じquadricの平方）。
正規化済みの全源の自由係数はn₀,n₃,n₄の三個で、失われたn₅は
実変換 `G=I+6n₅N²` の逆で一parameterとして戻る。
従って固定div P=Nの全対象支持は次元高々四。
`div P=tN,t≠0` も許すcone Tはtensorのscaleを一つ加えて次元高々五である。

GL₄内でNと可換な
`G=a₀I+a₁N+a₂N²+a₃N³`, a₀≠0、という四次元群KはTを保つ。
PGL₅×P(T)の像へのfiberには `(g,P)→(gk⁻¹,kP)` という自由K作用があるため、

`dim 全PGL像 ≤24+(5−1)−4=24`。

強い局所唯一性や全成分認識は不要である。
これはN²を使ったnormalizationを形式因子変換だけで済ませず、
旧証明書の全六括弧一致を入力にする全source次元論である。

## 3. 型5：quadratic不変核とx²係数が矛盾する

斉次場は `N=y∂x+z∂y+h∂z`。
その全不変二次式は正確に

`span{h²,hw,w²,2hy−z²}`。

nilpotent場で安定な各既約因子は固有値零なので、qᵢもこの核に入る。
従ってq₁q₂はxに依存しない。

一方、`HYPERPLANE_QUARTIC_NIL31_LONG_SUPPORT.md` の全Jacobi支持で
n₁₄=a≠0の二完全源の実Pfのx²係数はいずれも3a/2≠0である。
従ってその全点は対象にならない。
a=0の全補集合でのPfは、同ノート第1節により
線形因子の平方を持つもの、quadric平方、最高quartic零のいずれかである。
いずれも相異なる二既約quadricと単純Hの型ではない。
この排除では、既約quartic源の一般点だけを調べていない。

## 4. 型8：最小nil証明に必要な仮定を正確に弱める

N=y∂xとする。`HYPERPLANE_QUARTIC_MINIMAL_NILPOTENT_SUPPORT.md` 第1節で
Fの既約性から最初に必要としたのは `F|y=0≠0` だけである。
二つの既約quadricの積もyで割れないので、この仮定を満たす。

同節の全modular/Casimir恒等式から `F₀=a₀f₀` は非零定数であり、
a₀,c=f₀は非零定数、`l=αy,m=βy` を得る。
α=β=0ならFはyだけの多項式となり、斉次化後に二既約quadricを持てない。
その他の全支持は実GL₂変換でl=y,m=0へ入り、同ノート第2節の完全PDE核1,Iから

`F=a₀c−I`,
`I=y³w+c y²/2+u y³/3+v y⁴/4`。

modular恒等式によるIの係数−1の強制もF≠0だけを使う。
a₀c≠0なのでこれはwについてprimitive一次式であり、絶対既約quarticである。
従って二quadric積とはならず、全型8の対象支持は空である。
既存独立PDE/minor・全source reviewerの算術を再利用し、
旧MINIMAL_M成分の認識そのものはこの排除に不要である。

## 5. 型9：非零u係数は可約なので全HQQは平面となる

`Z=y∂x+∂z`、斉次化した `u=hx−yz` を使う。
全不変二次式は

`q=a u+B(h,y,w)`。

`HYPERPLANE_QUARTIC_NIL211_SHORT_TRANSLATION.md` 第1節の全potential表示は
quartic既約性を使わず、全Jacobi上で

`F=−u c−φ+k`、
`c=n₁₈h²+n₁₉hy+n₂₀hw−n₂₅y²/3`、`n₂₀n₂₅=0`。

n₂₅=0ならcはhで割れ、n₂₀=0ならcは(h,y)のbinary quadraticなので、
全Jacobi支持のcは常にrank二以下である。

c≠0と仮定する。Fはuについて一次で、二quadric因子のu係数はそれぞれ定数だから、
ちょうど一方の因子だけがuを含む。
uを含まないもう一方は、u係数cに比例する。
これはrank二以下で可約となり、qᵢの既約性に矛盾する。
従って必ずc=0。実Dは(h,y,w)だけに依存する平面型であり、
既済の全七conic配置・全磁気sourceの次元論により全像次元≤24。
個別の平面源の既知成分帰属をこの議論では仮定しない。

## 6. 型10：quadratic×quadraticの全枝は三parameter族

`HYPERPLANE_QUARTIC_TRANSLATED_MINIMAL_NORMAL_FORMS.md` 第1–3節の
全25係数potential還元は、quartic既約性を課す前に

`F=δR²/2+H`, `R=x−y²/2`, `δ≠0`,
`Vbase(H)=−H`, `tr Vbase=−1/2`, `deg H≤3`

へ帰着する（L除去は実affine変更）。
非零Hの全weightは同節の三型に尽きる。
H=bz又はbz³、及びH=z²(a+bw),b≠0では既済のGauss証明からFは絶対既約quartic。
H=0なら同じquadricの平方なので除く。
従ってHQQが残るのは `H=az²,a≠0` だけである。
この枝では零weight方向のtranslation κ∂wは任意に保つ。

canonical全族の自由係数は `δ,a∈C*`, `κ∈C` の三個であり、全括弧は

`Pxy=R/2`, `Pxz=−yz/2`, `Pxw=−2az+κy`,
`Pyz=−z/2`, `Pyw=κ`, `Pzw=δR`。

全Pfは `δR²/2+az²`。
`z'=g z,w'=d w+f`, gd≠0、という三次元affine群Kがこの全族を保ち、
`δ'=gdδ`, `a'=da/g`, `κ'=dκ`。
PGL₅×(三parameter族)への自由fiber作用から、全射影像次元は

`24+3−3=24`。

tensor全係数を変換したcanonical族の表示を使うため、
Pfの因子型だけの一致から成分を推測する問題はない。

## 7. 検証の区別と残る義務

第1表の全線形源・Jacobi・Pfは旧独立90係数reviewを入力とする。
regular及びnil31-longの**全非線形支持還元**はそれぞれの旧有限算術証明書と
記述証明を入力とする。
新しい独立target `code/classification_two_quadrics_review.py --nilpotent-factor-reductions`
は、型5/9/10の完全不変二次核、型1の全HQQ stratumの元Jacobi/Pf恒等式、
型3の無条件z因子、型9の全28parameter potential変換・全Jacobi・実係数n₂₀n₂₅、
型10の全canonical族と三次元群作用をFractionだけで再構成してPASS（0.023秒）。
証明書は `CLASSIFICATION_TWO_QUADRICS_PURE_NILPOTENT_FACTOR_REVIEW_20261107.json`、
実行前gateは `jobs/two_quadrics_pure_nilpotent_factor_review.json`。
型8の全PDE核及び全sourceは既存独立reviewを再利用する。
旧理想所属の全証明を新たに独立再構成したと混同しない。

本稿のG用次元論と、各源のE成分認識は別である。
型4の独立 `CLASSIFICATION_HYPERPLANE_QUARTIC_NIL31_SHORT_REVIEW_V2_20261107.json` と、
型7の独立 `CLASSIFICATION_HYPERPLANE_QUARTIC_NIL22_LONG_SUPPORT_REVIEW_20261107.json`
により、全十型について必要なG次元上界が揃った。E側の未帰属源は別に残る。
