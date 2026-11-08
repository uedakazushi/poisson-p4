# K4a actual24 親＝T37：可逆 PGL・Euler 全恒等式の独立 review

2026-11-08、PRIMARY `E_reduced_nonnormal`。
**別作者の QQ[θ]/(θ²+6) Fraction 算術と whole-source の幾何論理は全て PASS。**
K4a の稠密 actual24 軌道とその親は既知 T37、元の全 K4a source も T37 に属する。
新成分の増分は零である。
他の weak wall、全 E/G 分類又は総成分数をこの一つの同一性から決めていない。

## 1. 真の完成 source と同じ actual target

source は `MPS_SMOOTHING_COMPLETION_20261103_v62/T37.json` の
p0=t0=t1=t2=t3=1 とした全 completed quadratic tensor。
七つの高次補正を一つも省かない。
その五座標 x0,…,x4 を X=(h,x,y,z,w) と書く。
この全族の T37 component への帰属は accepted MPS full-family theorem を用いる。

target は accepted K4a actual polynomial curve の e=1、元 graph parameter t=174。
Y=(c,a,b,u,v) で、全十括弧を
`CLASSIFICATION_K4A_ORDER_TWO_ZERO_LOCUS_20261108.json` から読み直した。
先行 order-two independent review はこれが全 graph の同じ evaluation と確認している。
新しい point、frame、normal、profile 又は rank を選び直していない。

独立実装 `classification_k4a_t37_actual_map_review.py` は、standard-library `Fraction` の対
(a,b)=a+bθ、θ²=−6 の加減乗除と五変数疎多項式を独自に実装した。
producer、current helper、SymPy、CAS を import・実行しない。
全入力 SHA を固定し、active finite job を先に gate PASS とした後、一回実行した。

## 2. 可逆写像、非零 scale と全 Euler 項

次の X→Y の matrix を literal に独立実装した。

```
M = [[ 2, -6,  6, 3, -3],
     [ 3,-18, 18, 6, -6],
     [ 0,-3θ,-3θ,2θ, 2θ],
     [ 0,-3θ,-3θ, θ,  θ],
     [11,-54, 54,18,-18]]
```

QQ(θ) の Gaussian elimination で M^{-1} と det M=2592 を再構成し、
M M^{-1}=M^{-1} M=I を全二十五 entry ずつ確認した。
λ=−θ/36 は非零、λ²=−1/216。
θ のどちらの complex root でも全証明が成立し、未処理の unit/zero 補集合は無い。

Y 座標で全一次 vector field を

```
W/θ = (u−7b/9, 2u−3b, 3c−13a/4−7v/4,
       3c+5a/4−3v/4, 7u−5b)
```

とする。元の全括弧から独立に matrix multiplication／polynomial substitution を行い、

`M P_T(X) Mᵀ − λ P_K(MX) − (MX)∧W(MX) = 0`

の十 skew entry ×全十五二次 monomial、**全百五十係数**を確認した。
さらに X=M^{-1}Y を別に代入し、Y 変数の

`M P_T(M^{-1}Y) Mᵀ = λ P_K(Y) + Y∧W(Y)`

も全百五十係数で再確認した。
forward の座標方向や scale を逆に読んだ等式ではない。
全 M^{-1}、M、W と λ は独立 receipt に保存した。

Euler 項 Y∧W は projective bivector では零。
従って Φ:[X]↦[MX] は実際の PGL5 同値を与え、
Φ_*[P_T]=[P_K] が射影 Poisson parameter 空間の恒等式として従う。
radial 項を捨てることが出来るのは、この全 Euler 項を係数で確認したためである。

## 3. 全 signed Euler Pfaffian の独立 readback

各元 P から

`F_P = Σ_i (−1)^i X_i (P_ab P_cd−P_ac P_bd+P_ad P_bc)`

を定義通り独立に構成した。
ここで a<b<c<d は i の補集合。
全百二十六 degree五 monomial に対して

`det(M) F_T(X) − λ² F_K(MX) = 0`

を確認した。従って F_K(MX)=−559872 F_T(X)。
normal symbol、部分 Pfaffian、因子型又は零集合だけの一致ではない。
全 F_T と F_K も独立 report に係数として保存した。

全三百五十 projective Jacobi の零は、先行の真 completed source と
accepted actual graph の正式独立 receipt を利用する。
今回の可逆 actual PGL・scale・Euler 恒等式はその構造を保つので、
新しい Jacobi operator や matrix を再計算する必要は無い。

独立算術は **0.43508秒、RSS12,160KiB**で全 PASS。
正式 receipt は `CLASSIFICATION_K4A_T37_ACTUAL_MAP_SEPARATE_REVIEW_20261108.json`。
実行時 checker の全 bytes と SHA を埋め込み、source 二つと exact field のみから再生出来る。

## 4. 一つの実同一性から全 K4a へ戻す矢印

accepted `K4A_COMPLETE_TORUS_SLICE_SEPARATE_REVIEW.md` は、
同じ canonical P0 に対する全 e の actual curve、
全 e≠0 が一つの actual24 orbit、
その閉包 C が P0 の唯一の二十四次元 component と証明している。
この theorem の D101／G48／IFT graph 又は rank を今回は再計算していない。

上の全 PGL 恒等式から P(1) は既知 closed PGL-stable component T37 に属する。
従ってその actual24 orbit 全体と閉包 C は T37 に入る。
両者は既約で dim C=dim T37=24 なので C=T37。
実曲線の閉性から P(0) も T37 に属する。

さらに `PLANE_CONCURRENT_M3_NONZERO_QUARTIC_SEPARATE_REVIEW.md` の既受理全 K4a 正常化は、
β=1 の全 magnetic 補集合を factor で排除し、β≠1 を全 actual shear で戻し、
全 maximal-pole recurrence の n=0/1 補集合を保持した上で、
n=2 の全 A非零、B、linear magnetic 項を実逆付き shear/scale で P0 の一軌道へ戻す。
この whole-source theorem と T37 の PGL 不変性により、**元全 K4a source が T37 に属する**。
chosen magnetic slice や標本一点から全 parameter へ外挿していない。

以前の order-two-zero による T36/T38 陰性は整合するが、今回の陽性証明の仮定ではない。
T37 generic source の diagonal torus 正常化、全零点分類、intrinsic frame 構成の一意性、
discovery appendix の全係数探索網羅性も成功証明の仮定ではない。
必要なのは二つの公開 actual tensor と一つの可逆 matrix の全恒等式だけである。

本 receipt は K4a の named-parent 認識を完了する。
他の concurrent m三 branches の既受理 E認識と組み合わせる統合は、その既存被覆定理の範囲で行える。
新成分計数、two corrected-triple weak walls の既知親、全研究完了は主張しない。
