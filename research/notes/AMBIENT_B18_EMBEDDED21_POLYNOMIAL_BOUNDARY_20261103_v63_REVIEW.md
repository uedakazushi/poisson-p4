# embedded21 多項式境界の独立監査

2026-11-03。対象は `AMBIENT_B18_EMBEDDED21_POLYNOMIAL_BOUNDARY_20261103_v63`。
既存の W≠0 の実際の交叉認識を入力として受理し、多項式境界への輸送を別に監査する。
新しい source jet、Gröbner 基底、原点の交叉環計算は行わない。

## 1. 独立な係数確認

変数を `(h,x,y,z,w)` とし

```
E(W,λ)=∂y∧(x²∂h+(hx−λh²)∂z+(h²+z²)∂w)+W x²∂z∧∂w
```

を直接作る。全25発散係数、350生 Jacobi 係数を QQ[W,λ] 上で検査する。
固定 N19 の基底を、二つの f 方向、一つの g 方向、16個の Hamiltonian
方向から新たに再構成し、正準 Euler 商を取って全150係数で照合する。
保存済みの19係数を掛けたものは E に等しい。

対角行列を直接代入すると

```
r^−1 diag(1,1,r,1,1)_*E(W,λ)=E(W/r,λ)
```

となる。点 `[0:1:0:0:0]` では Euler 方向は x であり、その行列の x 行・列は
零である。残りの4次小行列式は W²。従ってこれは実際の射影 Poisson rank4
の証明である。W=0 では ∂y で分解し、この点の rank は2なので一般 rank2。

旧自由四係数 potential を直接作り、保存済みの frame を逆輸送する。
k0=1, k1=λ, k2=0, k5=1/3 で E(1,λ) に等しい。author の正規化関数は呼ばず、
不変の v02 JSON の free61、eliminated55、Euler 商から元の150係数を復元する。
その正準像がこの potential の bivector に一致し、perspective と vertical の
全 pivot は零、value は固定値となることも確認する。

旧49標的での値は、共通10座標が

```
(0,−k1,0,0,k0,2k2,0,3k5,0,0)
```

で残る35座標が零、first jet は固定である。これは旧 embedded colon prime
の自由四座標の族そのものである。単に同型な別の模型との照合ではない。

P0=E(0,0) の19個の GL 方向は係数 Lie 作用

```
δPij=δia Pbj+δja Pib−x_b∂aPij
```

から独立に作る。∂W, ∂λ の二列と合わせ、全150×21成分を照合する。
保存された21行による小行列式は QQ 上で32である。

## 2. associated support と全多項式境界

入力となる独立受理済み W1 認識は、実際の D=CA∩CB の完備局所環に
極小 prime と一つの正則21次元 embedded prime があることを証明している。
忠実平坦な完備化の associated prime の公式により、その収縮は実際の
代数的 associated prime である。excellence と完備化の次元公式により
その閉包 S21 の次元は21。被約 N の極小 prime とは次元で異なる。

G=GL5×tensor-scalar は連結で D のイデアルを保つ。G は同じ次元21の
associated support の有限集合を各々保つ。具体的に、一つを V とすると
G×V の像の閉包は既約で、有限個の21次元 associated support の和に入る。
従って一つに入り、恒等元の像 V を含むためそれは V 自身である。
高次元の極小 support に V が入ることは、この同次元の議論を妨げない。

上の自由四変数族の W1 での形式近傍は、実際の embedded prime に入る。
S21 の方程式をこの多項式族へ引き戻すと、局所完備化で零になる。
多項式環からその一点での完備化への単射性により、全引き戻しが零である。
従って全四変数族が S21 に入る。固定 frame と G 不変性から全 E(1,λ) が、
さらに対角輸送から W≠0 の全 E(W,λ) が入る。閉性により W=0 も入り、
特に P0∈S21 である。新しい associated prime を別に数えてはいない。

## 3. 原点で受理する形式的結論

GL19 局所切断×A² の21変数写像は P0 で小行列式32を持つ。
形式的逆写像定理により、選んだ21標的座標は source の形式座標であり、
この写像は正則21次元の形式的閉像を与える。像は S21 に入り同次元なので、
S21 の完備化の一つの極小形式枝である。共通 G19 を除けば正則2次元の
形式的像である。この主張は S21 の全 germ の正則性や解析的既約性を含まない。

同じ rank21 と既約な S21 の次元21は、21変数の支配的有理写像を与える。
ここで「有理 parametrization」はこの意味に限る。双有理性、S21 自身の
有理性、global proper resolution は導かない。

P_a、Q_b による embedded20 support とは次元で異なる。これにより原点の
nil support が曲線の横断部分だけでは尽くされないことが分かる。ただし
原点の D 全体、nilmodule、深さ、associated prime の網羅、新しい Poisson
成分については結論しない。

## 4. 再現

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_B18_embedded21_polynomial_boundary_review_20261103_v63.py
```

測定時間、最大 RSS、入力 hash と独立証明書は同名 REVIEW report folder に保存する。
legacy/v02 の三つの JSON は保存済み `sources/v02.zip` の展開物であり、
両方の hash を記録する。author source を import せず、既存レビューからは
小さな Fraction 行列式ルーチンだけを用いる。
