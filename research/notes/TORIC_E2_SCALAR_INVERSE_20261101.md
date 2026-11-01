# トーリック三次被覆の2-torsion座標：逆写像と trace

2026-11-01。`TORIC_JACOBIAN_TWO_TORSION_20261101.md` の代数同型に対し、
明示的な逆写像と普遍的 trace 恒等式を与える別の証明書である。
トーリック成分の一般三次被覆を扱い、全 S61 ファイバーの尽くしは主張しない。

## 1. 三次共変量 D

primitive 閉一次二形式を β、`T_a=(β_ij,a)` とする。添字は0から3。
標準 Poisson 行列を J とし、対称行列

\[
 M_{ab}=\sum T_a{}_{ij}T_b{}_{kl}J_{ik}J_{jl},\qquad
 N_{ab}=\sum T_j{}_{ai}T_l{}_{bk}J_{ik}J_{jl}
\]

を定める。各和は i,j,k,l の全範囲を取る。
対称三次形式 D を

\[
 D_\beta(y)=\frac92\beta_y(y,JMy)-36\beta_y(y,JNy)         \tag{1}
\]

と置く。その係数は β の三次式である。
対称テンソル成分 D_ijk と多項式係数を区別し、例えば y_i²y_j の係数は3D_iijとする。
交代的な不変縮約の順序を

\[
 \langle D,c\rangle=\sum D_{ijk}c_{abl}J_{ia}J_{jb}J_{kl} \tag{2}
\]

と固定する。

(1) は自然な tensor 縮約なので Sp4 同変である。さらに全16変数を保ち、
Lie 代数の全10方向について20成分、合計200本の恒等式を厳密に確認した。
`code/toric_degree3_scalar_20261101.py` はこの計算と、二つの既知の有理標本における正規化を保存する。

## 2. 普遍的縮約と稠密断面

U=C5 は既存の五次共変量、V=C9 は結合則の線形部分 B による

\[
 B(V)=-H(U,U)/3-27c_4C                                  \tag{3}
\]

の一意な一般有理解とする。V の存在・全結合則との整合性は主証明書を用いる。
B の一般階数20は本コードでも別に再構成した有理点で確認した。
(3) の自然性と単射性から V の Sp4 同変性が従う。

本証明書は次の両式を独立に確認した。

\[
 \langle D,U\rangle=54c_4,\qquad \langle D,V\rangle=162c_6.\tag{4}
\]

これは標本一致からの推測ではない。重み座標 b の六次元斉次線形断面

\[
 b=s_0(1,2,\ldots,16)+s_1e_{11}+s_2e_{12}+s_3e_{13}+s_4e_{14}+s_5e_{15}
\]

へ D,U を独立に制限し、QQ[s₀,…,s₅] で確認した。c4,c6 は同断面で Fisher の公式を直接計算した
427項・2090項の多項式を使う。
断面と10個の Sp4 接方向の16次行列式が非零であることも独立に再計算した。
従って断面の Sp4 軌道和は W16 に稠密である。

両辺は Sp4 不変である。右辺については canonical Pfaffian curve の変数変換・合同変換の
Fisher character が Sp4 上で1となる。V の有理定義域を縮小すれば、稠密性により (4) は
一般 W16 上の恒等式となる。

## 3. 三次代数の逆写像

K=QQ(W16) とし、判別式を可逆化する。対象は canonical genus-one curve の
Jacobian の非零2-torsion代数

\[
 A_E=K[X]/(f),\qquad f=X^3-27c_4X-54c_6
\]

である。Cβ 自身には一般に選ばれた原点がないため、`Jac(C_beta)[2]` と記す。
主証明書の三つの96成分恒等式は

\[
 c_i=(U_iX+V_i)/f'(X),\qquad f'(X)=3X^2-27c_4             \tag{5}
\]

が全結合則を満たすことを示す。(4) と `Xf'−54c4X−162c6=3f` により、逆写像は

\[
 \boxed{X=\langle D_\beta,c\rangle}.                     \tag{6}
\]

既存の全96式の一般商環は次元3なので、(5) が与えるその商環から A_E への写像は
(6) により全射、同次元により同型である。この箇所では「全96式の一般商環が次元3」という
以前の定理を用い、特定の三解だけから核を零としていない。

分母の可逆性は a=27c4,b=54c6 と置いた形式的 Bézout 恒等式

\[
 (3X^2-a)(6aX^2-9bX-4a^2)
 -(18aX-27b)(X^3-aX-b)=4a^3-27b^2
\]

から従う。右辺は c4³−c6² の非零定数倍である。

## 4. 普遍的平均の恒等式

形式三次代数の乗法行列を用い、a,b を独立変数のまま

\[
 \operatorname{Tr}(1/f')=\operatorname{Tr}(X/f')=0
\]

を確認した。(5) により全20座標の trace が零となる。
普遍的アフィン復元 `q=bar q+N c` と合わせると

\[
 \boxed{\operatorname{Tr}(q_T)=3\bar q=2Q_{51}+Q_{52}}.   \tag{7}
\]

従って以前の標本での三トーリック点の和は、一般一次部分上の普遍的恒等式へ強まる。
この証明は余次元1における有限性や、C9 の全 W16 上の多項式延長を必要としない。

## 5. 再現と依存

```sh
.venv/bin/python code/toric_degree3_scalar_20261101.py \
  --output reports/new-run/toric-degree3-scalar.json
.venv/bin/python code/toric_e2_scalar_inverse_20261101.py \
  --output reports/new-run/toric-e2-scalar-inverse.json
```

証明書は `TORIC_DEGREE3_SCALAR_20261101.json` と `TORIC_E2_SCALAR_INVERSE_20261101.json`。
後者は稠密断面上の縮約、断面の優越性、B の一般単射性、形式的逆元と trace を確認する。
全96結合則の residue 恒等式は主証明書への明示的な依存である。
平面外の未知の S61 解や判別式零上の特殊ファイバーは、この結果では排除・分類しない。
