# 二つの Y wall：全 higher source に既約五次因子はない

2026-11-08。PRIMARY `E_reduced_nonnormal`。producer 全源証明。独立確認は別工程。
新 M1、点、周辺 tangent profile、minor 選別、jet は使用しない。

## 1. どの未処理全源を閉じるか

corrected V2 の全 first-connection 分類における、直接 G 次元上界だけで閉じていた
`E1_LOG_GENERIC_Y, τ=−1` と `E2_LOG_GENERIC_Y, τ=−4/3` を扱う。
両方とも保存済み full350×30 weight 一 operator の全 kernel は九次元、
full350×9 weight 二 operator は injective である。
**injective は P2=0 を意味しない**。全九つの inhomogeneous P2 を必ず戻す。

保存済み完全核 K1 の全九係数を p0,…,p8、全 P2 を r0,…,r8 とし、
`Π=P0+P1(K1*p)+P2(r)` とする。P1 は normal-weight 一の全三十 source、
P2 は normal-weight 二の全九 source からそのまま作る。
K1 による leading actual D3 の直線条件も元のままである。
全 normalizer 逆は corrected V2 first-symbol proof のものを使い、
canonical higher source の任意の一点へ縮めていない。

## 2. 全 source の literal Euler/Ore 変換

座標(a,b,c,z,w)。両 complete K1 で、P1_ab の c*z,c*w coefficient に対応する
行2,11は全九列で零である。他の non-c entry は source 定義から c 非依存、
全 P2 も c 非依存である。また Π_cj は c に関して高々一次である。

`R_c=0, R_j=∂c Π_cj (j≠c)` とし、実 Euler quotient の同値代表

`Π'=Π−E∧R`

を取る。R は c 非依存の一次 field であり、Π' は **全 source で c 非依存** となる。
全350 projective Jacobi polynomial の Euler correction 前後の恒等一致を確認した。
従って `Π'=Q(a,b,z,w)+∂c∧W(a,b,z,w)` である。
元 `E∧Jac(Π')` の c 一次係数は符号だけを除き全 strict Jacobi(Q) そのものであり、
actual Poisson source では strict QJac=0 が必要である。変換で新しい stabilization 仮定は置かない。
全 Q/W は report に保存した。

P0 だけについては R_w=5z/3（第一）、z/3（第二）であり、
Euler correction 後の core はそれぞれ

```
第一: Qaw=−4az, Qbw=−bz, Qzw=0,
第二: Qaw=−2az, Qbw=bz, Qzw=3z².
```

しかし証明にはこの P0 だけでなく、全九 P1 と全九 P2 を使用する。

## 3. 全 inhomogeneous P2 の正確な消去

保存済み weight 二 matrix を B、全 P1 self-Jacobi を C(p) とすると、
weight 二の元方程式は `B*r+C(p)=0` である。
保存済みの同じ rank 九 minor（固定 QQ 非零定数）の九行から一意の
`r=r_part(p)` を解き、残350 compatibility と全 original350 Jacobi を全て残す。
未知 source coefficient の逆元、generic specialization、homogeneous P2 だけへの制限はない。

full original350 の monomial order は helper の五quad、degree四 exponent 降順である。
P2 を上の一意の quadratic polynomial で置換すると、
**元の row54 と row279** は次の値になる。

```
第一: J54=J279=2*p8²,
第二: J54=J279=2*p8²/3.
```

従って C 上の reduced Poisson support では必ず p8=0。
非被約 scheme の nilpotent 指数を求める必要はない。

## 4. actual divisor の全因子恒等式

全 original actual divisor を `D=Σ_i (−1)^i x_i Pf(Π_omit_i)` とする。
P2 の一意式を代入し、必要条件 p8=0 を戻すと、正確に

```
第一:
D=−z³*(-4a+p5*w+p7*z)*(-b+p3*z+p6*w)*(p0−4p1−p7)/4,

第二:
D=z³*(p0−p2)*(-5a+p1*z+p5*w)²/5.
```

全 source parameter の polynomial identity であり、先頭係数が零の場合も含む。
いずれも D=0 又は proper linear factor を持つ reducible/nonreduced divisor である。
従ってこの二つの全 canonical higher source は rank 四・被約既約五次 D を持たない。
全 first-symbol normalizer の真の逆を戻しても既約性は保たれるので、
対応する二つの全 Y-wall stratum は全点で空である。
G 次元上界だけの主張を、この二 stratum について E 全源空性へ強めた。

第一の full D は p8=0 を使う前から factor `b−p3*z−p6*w` を持つが、
両源に共通な original square-row と上の完全な積式を主証明とした。
他の constant-root/moving-root/triple first-symbol branches は別証明のままであり、
この二源だけから全 E または全分類の完了を主張しない。

## 5. 証明書

事前 gate PASS、180秒/3GiB。新計算は保存済み九核の Euler/Q/W 再構成と、
同じ固定九 minor の inhomogeneous elimination、元全350と actual D の恒等比較だけ。
5.685秒、最大RSS246540KiB。

`CLASSIFICATION_PLANE_TRIPLE_Y_WALLS_ORE_20261108.json` に全 P1/P2、R、Q/W、
固定九 minor、一意 r_part、全350 compatibility、全350 original Jacobi、full D を保存した。
旧五十二 FunctionDef AST、58c43 source、当時の全41 binding契約は preappend archive に保持した。
昔の G/E report、producer source、manifest、current binding は書き換えていない。
