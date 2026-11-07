# contact4_secant_at_contactの全六磁気源はβ=0へ移る

2026-11-07、課題E_reduced_nonnormal。
`TWO_QUADRICS_PLANE_CONIC_CONFIGURATIONS.md` の七free配置の一つを、
**全六磁気係数で**既存β=0支持へ認識する。
この源の二次・三次法線支持は非孤立だったが、そこから実変形を推測せず、
全源の正則な座標変更を直接構成する。

## 1. 完全源と奇数次数項

base座標をh,x,y、fiberをz,wとすると、二conicは

`q₁=xy−h²`, `q₂=xy−h²+x²`。

完全Saito anchorは

`U=−x²∂x/2+(h²−xy/2)∂y`,

`V=(10h²/9−10x²/9−5xy/3)∂x+(10xy/9+5y²/9)∂y`。

`[U,V]=aU+bV`, `a=10x/9+5y/3`, `b=−x/2` であり、全源は

`P=U∧∂z+V∧∂w+(-az-bw+R)∂z∧∂w`、

`R=m₀h²+m₁hx+m₂hy+m₃x²+m₄xy+m₅y²`。

全mᵢは任意であり、実Pfaffianは `10q₁q₂/9`。
h=1での一次項は正確に `(m₁x+m₂y)∂z∧∂w` だけである。
U,Vはaffine baseについて偶数次数、a,bはx,yの一次式である。

## 2. 全parameterに正則なfiber translation

`z'=z+αh`, `w'=w+βh` と置くと磁気項は `R'=R+h(αa+βb)` へ変わる。
α,βからhx,hy係数への写像の行列式は5/6で非零。
従って

`α=−3m₂/5`, `β=2m₁−4m₂/3`

とすれば全一次項が消える。
全六parameterに多項式的で、parameterの反転や除外はない。
この操作はconic二式を変えず、元のh-Casimir strict Poisson全源を保つ。

translation後の定数bivectorCは

`C^{xw}=10/9`, `C^{yz}=1`, `C^{zw}=m₀`、他は零。

Pf C=10/9なので、全mᵢについて非退化である。
さらに一様なDarboux線形変更

`(X,Y,Z,W)=(x,9w'/10,y,z'−9m₀x/10)`

はCを標準 `J=∂X∧∂Y+∂Z∧∂W` へ送る。
全座標変更のdeterminantは9/10であり、どのparameter値でも可逆。
線形変更は消えた一次項を再導入しないので、**全源が正確に**

`P_new=h²J+Q(X,Y,Z,W)`, `[J,Q]=0`, `[Q,Q]=0`

となる。β=0の接方向だけを確認したのではなく、全tensorをこの形へ移した。

## 3. 全支持への結論と独立算術

`BETA_ZERO_COVERAGE.md` の既済定理は、この全Q45支持を既知52成分の有限閉和へ入れる。
従って本源の全六磁気係数とPGL飽和も、その既知閉和へ入る。
特定一点でのcomponent membership、局所唯一性、又は法線pencilの存在は不要である。
どの既知成分番号へ入るかの細分類を、この被覆に追加仮定としない。

再利用producerの `--plane-even` modeは全座標行列・全六変換後括弧・C=J・
全一次係数零・全Jacobi・Q45互換条件を文字通り確認した（0.124秒）。
独立 `classification_nil22_base_review.py --plane-even` はproducer/CASをimportせず、
逆affine置換と六括弧のpushforwardを別のFraction疎多項式で構成し、
全parameterで同じ恒等式を確認した（0.004秒）。
専用gate通過後の証明書は
`CLASSIFICATION_TWO_QUADRICS_PLANE_EVEN_20261107.json` と同名 `EVEN_REVIEW`。
既存β=0の尽くし自体と新しい全源座標恒等式の確認は区別する。

残る未帰属平面源はcontact31_tangent_at_simple_intersectionとcontact31_chordの二つ。
両者の既知T39/T40全graphへのwitness検査は陰性だった。
その陰性から新成分性や既知成分への非包含を一般には推論せず、全source構造を別に調べる。
