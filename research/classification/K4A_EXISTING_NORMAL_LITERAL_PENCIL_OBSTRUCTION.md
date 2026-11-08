# K4a の保存済み唯一 normal character：一つの literal pencil の非閉鎖

2026-11-08、PRIMARY E。新しい point/operator/rank/minor/normal Jacobi map は計算しない。
独立受理済み Dkernel49 と actual M51 の、同じ48行48列証明書を読むだけの有限 job を実行した。
新しい miniversal/jet 延長を行わない。

## 固定行列からの正確な quotient 読み戻し

`CLASSIFICATION_K4A_CANONICAL_FULL_AMBIENT_20261108.json` の全 kernel は49列、
全 actual Euler/GL5/scalar像は48次元。
torus character -1 の既 kernel 列を順序通り読み、同じ既 M48 block で解いて全150列を比較した。
最初に非groupとなるのは zero-based **kernel列23 / 自由係数列112**。
残差、M48係数、全51 actual M列を殺す左汎関数を
`CLASSIFICATION_K4A_EXISTING_MATRIX_NORMAL_READBACK_20261108.json` に保存した。
その汎関数は当該 N に **8/87** を与えるので、非group性は厳密である。
実測0.059秒。inline 実装全文は finite job と報告に保存し、凍結runnerを変更しない。

座標 `(c,a,b,u,v)`、当該既 cocycle は全体として

```
Nca=-12bc/29, Ncb=-70ac/29,
Ncu=37a²/174-b²/29, Ncv=2ab/29-6cu/29,
Nab=-3a²/2+9b²/29, Nau=6bu/29, Nbu=au,
Nav=Nbv=Nuv=0.                                               (1)
```

全項が stabilizer torus の character -1 を持つ。
D*N=0 は既全49 kernel certificate に含まれるので、[P,N] の新operatorは作らない。

## 一つの N の自己 Jacobi を手で計算する

規約 `J_ijk=Σ_l(N_il ∂l N_jk+cyclic)` とする。
(1) の七式だけを使い三行を手で展開すると

```
Jcab= -216abc/841,
Jcav= 3a³/29+6ab²/841+36bcu/841,
Jcbv= 53a²b/841+18b³/841+6acu/29,
Jabv=0.
```

従って元の projective 条件、quadruple `(c,a,b,v)` の係数は

```
(E∧J)c a b v = c*Jabv-a*Jcbv+b*Jcav-v*Jcab
             = [34a³b-12ab³-174a²cu+36b²cu+216abcv]/841.      (2)
```

特に a³b の係数は34/841≠0。
これは strict Jacobi の非零だけでなく、Euler 商後の元projective Jacobiの非零である。
従って保存されたこの一つの N に対し `P+eN` は actual Poisson pencil でない。
同じ N からの高次補正、別 cocycle の列、別pointへ進めない。

## 射程

判明したのはこの literal pencil の不適合だけである。
K4aのambient reduced support、他の actual nonlinear family、既知24親への帰属、
又は新成分性を否定も証明もしない。
全D101・orbit23・tangent24という既 limitation はそのまま保持する。
次は既知24親の actual 構成又は別の source 幾何へ戻る。
