# E3の実軌道閉包：射影nilpotent接対称性による全境界の被覆

2026-11-07、課題G。主軌道の延長空性は既済であり、本ノートは実際の軌道閉包を扱う。
外部論文の積コンパクト化や有限個の対角1PSを全境界の被覆とは仮定しない。
**結論：E3成分閉包の全最高部について、全unimodular低次延長のprimitive β像は一般底を優越しない。**
証明は、射影的に閉じたnilpotent接対称性incidence、その全四Jordan型、
全80三次一形式係数の完全線形核、既済の閉対数源定理を使う。

## 1. 最高部と三次一形式の線形同定

V=C⁴、EをEuler場、固定体積形式をvolとする。
閉じた係数次数二の二形式γに対し

`θ=i_Eγ/4`, `dθ=γ`, `i_Eθ=0`                    (1)

である。逆に、Euler-horizontalな三次一形式θからdθによりγへ戻れる。
したがってこれは有限次元ベクトル空間の線形同型であり、非零点と射影閉包を保つ。
体積双対なbivectorをQとする。
`Pf Q=0` はγの可分解性を意味し、(1)はdegree2葉層の可積分条件へ対応する。
本ノートの全境界議論はこの実表現内で行う。

## 2. 非零nilpotent線形場は全極限で保持できる

`P(N)={ [Y]∈P(End V) : Y⁴=0 }` と置く。
これは射影nilpotent円錐である。特に全Yはtrace零であり、Y=0は射影化で除かれる。
θの射影空間との積内で

`I={ ([θ],[Y]) : i_Eθ=0, i_Yθ=0, Lie_Yθ=0, Y⁴=0 }`        (2)

を定義する。各条件は斉次多項式条件である。
Iのθ空間への射影像はprojectiveなので閉じる。
Yが退化してrankを落としても、[Y]はP(N)内で極限を持ち、非零性は失われない。

主軌道の最高部は

`Q₀=(5x²,5xy,5xz,y²+3xz,yz+7xw,7yw-3z²)`

であり、`Y=x∂y+y∂z+z∂w`、
`X₀=-5x∂x-y∂y+3z∂z+7w∂w` とすると `Q₀=Y∧X₀`。
`[X₀,Y]=-4Y`, `tr X₀=4`, `tr Y=0` なので
`div Q₀=0`, `Lie_Y Q₀=0`, `Y∧Q₀=0`。
trace零のYはvolを保つから、対応するγ₀についても
`Lie_Yγ₀=0`, `i_Yγ₀=0`。
YはEと可換なので、θ₀も(2)を満たす。
座標変更後もこの性質は保たれ、射影像の閉性から

`全E3軌道閉包 ⊂ pr_θ(I)`                         (3)

を得る。これはBorel軌道の任意のunipotent発散を含む実際の全極限についての包含である。
associated bundleや商空間のコンパクト化を実軌道閉包と同一視していない。

Cartan公式により、(2)は `i_Ydθ=0` も与える。
逆に `i_Yγ=Lie_Yγ=0` から(1)を使って(2)を得る。
よってbivector表示での正しい条件は `Y∧Q=0`, `Lie_Y Q=0` であり、
traceを落とした誤った双対化は使っていない。

## 3. 全四Jordan型の全80係数

非零nilpotent4×4行列のJordan型は(4),(3,1),(2,2),(2,1,1)の四つだけである。
各型で三次一形式の全80係数を取り、
`i_Eθ=0`、`i_Yθ=0`、`Lie_Yθ=0` の全係数を並べる。
これは150×80の有理数行列であり、いかなる係数sliceや可積分性条件も追加しない。
実際の境界を含む、より大きい線形源を分類している。

| Yの型 | 代表Y | 全行列rank | 完全核次元 |
|---|---|---:|---:|
| (3,1) | x∂y+y∂z | 73 | 7 |
| (2,2) | x∂z+y∂w | 73 | 7 |
| (4) | x∂y+y∂z+z∂w | 77 | 3 |
| (2,1,1) | x∂w | 65 | 15 |

以下は全kernelの数学的記述である。
producer reportは全行列、記述したbasisの全係数、相補的な非零rank minorを保持する。
独立checkerはこの数学的basis自体を別構成し、全行列と二つのminorを再構成する。
従ってbasisの完全性は標本やCASの未検証nullspace宣言に依存しない。

## 4. rank二の二型は全て閉L112へ入る

(3,1)では `r=x`, `s=w`, `P=y²-2xz`、
(2,2)では `r=x`, `s=y`, `P=xw-yz` と置く。
いずれもr,s,PはY不変であり、次数は1,1,2である。
第3節の完全核は次の七次元線形空間に等しい。

`θ=P(l₁ dr+l₂ ds) - (r l₁+s l₂)dP/2 + H₂(s dr-r ds)`,       (4)

ここでl₁,l₂はr,sの任意の一次式、H₂は任意のbinary quadraticである。
四つのl係数と三つのH係数を全て残している。
各式が(2)を満たすことは直接計算で分かり、相補minorはこれ以外の核がないことを証明する。

**補題。** r,sが独立、Pが上記のいずれかであるとき、(4)の全線形空間は
`overline L(1,1,2)` に含まれる。

**証明。** `C₂=-(r l₁+s l₂)/2` が相異なる二根を持つ条件は稠密開である。
r,sの可逆線形変更によって `C₂=krs`, k≠0とできる。
radial条件から `l₁=b s`, `l₂=c r`, `b+c+2k=0`。
従って

`θ=P(b s dr+c r ds)+krs dP+H₂(s dr-r ds)`。                (5)

binary quadratic `F₂` に対する線形作用素を

`F₂(b s dr+c r ds)+krs dF₂ = T(F₂)(s dr-r ds)`              (6)

と定義する。radial条件により右辺の形を持つ。
monomial `r^i s^(2-i)`, i=0,1,2への固有値は

`b`, `b+k`, `b+2k`。

これらが全て非零である稠密開では、`T(F₂)=H₂` を一意に解ける。
`q=P+F₂` とすると(5)は

`θ=rsq ( b dr/r + c ds/s + k dq/q )`                       (7)

となる。qが既約である条件も非空開である
（例えばH₂=0なら上記Pはrank3又は4の既約二次）。
従って(4)の非空稠密開は実L112源である。
七次元線形空間は既約であり、L112の成分閉包は閉なので、全(4)がその閉包に入る。
三固有値の零点、重根C₂、C₂=0、可約qを含む**全補集合**もこの包含が保持する。□

ここでTの可逆開のみを全低次延長の被覆と読み替えていない。
全最高部が既済の**閉L112**へ入ることを証明し、
`UNIMODULAR_LOGARITHMIC_POLE_BOUNDARY_REVIEW.md` により既にその全低次延長を処理している。
有限係数計算でも(6)の三恒等式を、自由なb,k,r,sに関する多項式恒等式として独立確認した。

## 5. regular nilpotent型はE3主軌道又は一次Casimir

Yをregular nilpotent `Nx=x∂y+y∂z+z∂w` とする。
第3節の完全核は、体積双対を戻すと

`Q=Y∧(a₀X₀+a₂N²x+a₃N³x)`                                (8)

の三次元空間である。
ここで `N²x=x∂z+y∂w`, `N³x=x∂w`。
全三次一形式80係数のrank77による証明なので、
pointwise分解を線形場のwedgeへ無根拠に取り替えていない。

a₀≠0なら、N²,N³はYと互いに可換であり、
`[N²x,X₀]=8N²x`, `[N³x,X₀]=12N³x`。
これらの線形場の指数によるGL4座標変更はYを固定し、
a₂,a₃を任意に消すことができる。高次交換子は零である。
従って非零scaleを除いてQはQ₀そのもの、すなわちE3主軌道へ入る。
主軌道上の全30低次係数の空性は
`UNIMODULAR_EXCEPTIONAL_SUPPORT.md` 第1節の独立受理済み定理である。

a₀=0なら、Y,N²x,N³xの全ての∂x係数が零なので
`Q^#(dx)=0`。
これは最高部の実際の一次Casimirであり、
`UNIMODULAR_LINEAR_HAMILTONIAN_SUPPORT.md` の全低次延長定理を適用する。
零Qも既済の零最高部に含める。

従って、実軌道の境界にregular型Yが残る場合はa₀=0へ入り、
境界から主軌道の空性だけを外挿することはない。

## 6. rank一型は三変数pullback

Y=x∂wとする。(2)とCartan公式から `i_Yγ=0`。
多項式環は整域なので `x i_(∂w)γ=0` は `i_(∂w)γ=0` を意味する。
さらにdγ=0から `Lie_(∂w)γ=0`。
従ってγはx,y,zだけから引き戻される二形式である。
`REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` 第16節の全pullback延長定理を適用する。
有限kernel15はこの全三変数水平一形式空間を確認する補助証明書であり、
この幾何的排除に特定のgeneric平面葉層の分類を使わない。

## 7. 実E3閉包と全最高部への結論

(3)と全四型により、最高部集合として

`overline E3 ⊂ E3_main ∪ overline L112 ∪ linear_Casimir ∪ pullback ∪ {0}`。       (9)

各右辺について**全低次係数を保持した**空性又はprimitive像の非優越定理がある。
従ってE3成分閉包の全延長は一般primitive底W16を優越しない。
この結論は外部のE3境界四成分記述に依存しない。

全最高部への外部入力は別に明示する。
[Pym, arXiv:1403.6444v1, Theorem3.1及び直後の段落](https://arxiv.org/html/1403.6444v1#S3)
を今回直接再読した。Theorem3.1は全可積分section Ω¹(4)の空間について
L1111、L112、R22、R13、S23、E3の六既約成分を述べる。
その直後は、CLNの元の余次元二特異集合の仮定に対し、
Loray–Pereira–Touzetの結果により、残りのPoissonに対応する葉層も
これらの成分閉包に入ることを明記している。
従って共通因子を持つθを黙って除外しない。
外部六成分定理とこのLPT拡張は外部入力であり、本有限checkerによる再証明ではない。

`REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` 第16節は既に、
R13閉包を閉一次Casimir条件、R22閉包をbasepointのないGrassmannian pencil写像、
S23閉包を閉pullback条件へ接続している。
L1111・L112の両閉包はquartic integrating-factor incidenceと全極型により処理済み。
最後のE3閉包は(9)で閉じる。
最高部γ=0とγ共通因子付きの場合も既済の全評価rank定理が保持する。
よって**A=c=0の全unimodular二次支持**のprimitive像はW16を優越しない。

他のモジュラー型は `REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` 第9・13–15節と
そこから参照する全源定理に従う。一部は全primitive像の非優越ではなく、全MPS帰属を与える。
したがってここで必要な最弱の結論は、**未知のW優越成分の一般元は重複超平面を持たない**
ということである。もし一般元が重複Hを持つなら、非優越枝を避けたその稠密開は
有限個の既知MPS閉成分の和に入る。未知成分の既約性と閉性から、その成分自体が
一つの既知成分へ含まれ、未知の全空間成分であることに反する。
追加の橋を示さずに、全重複H支持のprimitive像次元を一律15以下と強めない。
Gの完了を記録する際は、C0の全134式/S61・有限perspectiveの橋、
既知五点から未知成分への橋、正規/被約/非被約退化因子の全分岐が
現行計画で全て接続済みかを統合担当が最終確認する。
本ノートだけから全成分の上界89や特殊底の全被覆を主張しない。

## 8. 再現・独立検証

既存producer/reviewerの `--nilpotent-incidence` modeを使い、旧主軌道証明書は保持した。

```sh
python3 scripts/classification.py gate research/classification/jobs/unimodular_exceptional_boundary.json
.venv/bin/python code/classification_unimodular_exceptional.py --nilpotent-incidence
python3 scripts/classification.py gate research/classification/jobs/unimodular_exceptional_boundary_review.json
python3 code/classification_unimodular_exceptional_review.py --nilpotent-incidence
```

producerは6.50秒で完了。最初の実行ではforcing恒等式のSympy構文上の比較を
展開後の零判定に直してから再gate・再実行した。数学的源・条件・rankは変更していない。
報告は `CLASSIFICATION_UNIMODULAR_EXCEPTIONAL_BOUNDARY_20261107.json`。
独立結果は同名の `BOUNDARY_REVIEW` report。
独立実装はproducer/CASをimportせず、全80列、全行列、数学的に記述した全basis、
相補minor、(6)の三forcing恒等式をFraction疎多項式で再構成する。
projective閉性、稠密L112包含、regular conjugation、全最高部の橋は上記の数学的証明であり、
算術reportだけをその代わりにしていない。
