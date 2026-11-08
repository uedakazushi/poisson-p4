# 三変数pullback最高部：全unimodular低次延長の点別被覆

2026-11-08、PRIMARY `E_multiple_hyperplane`。producer純証明、独立/root査読前。
新arithmetical run、operator、minor、点、jetは使わない。

## 1. 元の全C/L分解を戻す

α=J+β+γは係数次数0/1/2の閉二形式、Pfαは非零constantとする。
非零constant vに対しi_vγ=0と仮定する。closednessからLie_vγ=0。
baseを(x,y,z)、v=∂wとしてγは全三変数pullbackである。

`REMAINING_DEGENERACY_DIVISOR_SUPPORT.md` §16の全低次分解を使う。
元の全一次項にはw dw∧c'が有り得るが、最高Pf条件からc'∧γ=0。
c'非零ならhighest dual Qがlinear Casimirを持ち、既に新しい全membership定理を適用出来る。
残c'零で

```
α=γ+β3+A+w Ω+dw∧(η+c), Ω=dη                    (1)
```

が全sourceである。β3はbaseの一次二形式、A/Ωはconstant二形式、ηはbaseのlinear一形式。
全Pf条件はb=η+cがbase全点で非零、b∧Ω=0を与える。
Ω非零なら、同節の全affine像階数の証明からη=ell*a、γ∧a=0となるnonzero constant aを得る。
従ってこの全補集合もhighest linear-Casimirの全membershipへ入る。
これら二つの枝をGの非優越だけで閉じたとは呼ばない。

## 2. Ω零ではfull tensorがconstant方向で不変

残Ω=c'零ではα全体がwに依存せず、bは閉affine一形式である。
従ってb=dH、Hはdegree≤2のpolynomial、dHはどこでも非零。
全inverse tensorもw非依存であり、固定orientationと非零constant scaleを除いて正確に

```
P=∂w∧B(x,y,z)+Jac(H),
div B=0, B(H)=k≠0, degB≤2.                       (2)
```

と書ける。Jacはbase順(x,y,z)のJacobian tensor。
(2)は元の全C/L/Qを保った表示で、PfP=k。
JacobiはLie_B Jac(H)=Jac(BH)−(divB)Jac(H)=0そのものである。

dHに零点が無いのでHess Hのrankは高々二、linear部はそのimage外にある。
baseのactual affine変換と非零constant scale/fiber補償により

```
H=z+q(x,y), rank q=0,1,2
```

の全三型へ移せる。Hのconstant差は元tensorを変えず、自由に選べる。
平方完成のimage内translationとradical-linear unitの逆を全て保持する。

## 3. rank q零／一

rank零ではH=z、B(H)=Bz=k。従ってP#dzは非零constant場であり、全CA membershipを適用する。

rank一ではH=z+x²/2と出来る。B2をBの二次最高部とすると、B(H)=kの三次係数は
x*B2(x)=0。整域上B2(x)=0である。
(2)の最高Qは∂w∧B2だけなのでQ#dx=0。
よって全linear-Hamiltonian最高部定理を適用して全sourceを有限既知閉和へ送る。
全linear/constant Bを残し、B2が零になる補集合も同じ定理に含む。

## 4. rank q二：全actual曲線で既NEW25へ

actual base frameでH=z+xyとする。Bに追加条件を付けず

```
B_e=B+e H∂z, P_e=∂w∧B_e+Jac(H)                  (3)
```

と定める。元座標でdegree≤2、eに関してpolynomial。
Hz=1から

```
div B_e=e, B_e(H)=k+eH,
Lie_Be Jac(H)=Jac(k+eH)−e Jac(H)=0.
```

従って(3)は全eで真のquadratic Poisson tensor、Pf(P_e)=k+eH。
e非零ではf_e=H+k/eを取れば、Jac(f_e)=Jac(H)、B_e(f_e)=e f_e。
すなわち非零constant modular場の**元全source**

```
P_e=∂w∧B_e+Jac(f_e), divB_e=e, B_e(f_e)=e f_e
```

である。f_eのquadratic Hessianはrank二、radical-linear係数は一で非零。
eでtensor尺度を補償し、baseのactual affine normalizationを戻すと
`CONSTANT_MODULAR_DEGREE_TWO_WHOLE_REDUCTION.md` のsource II、f=u+zwへ入る。
全C/Lはその必要十分な全A8源に含まれ、selected subfamilyだけを使わない。

`CONSTANT_MODULAR_SOURCE_TWO_NEW25_WHOLE_COVERAGE.md` と正式独立SEPARATE_REVIEWにより、
source IIの全parameter/frame/zero境界は既positive-weight NEW25に入る。
従って(3)の全e非零点はその閉成分に属し、e零の任意の元Pも同じ成分へ入る。
分母1/eをe零へ代入せず、元座標でpolynomialな(3)の閉性を使う。
k非零により曲線は射影点として消えない。
実因子は非零scalar倍h³*(h z+xy+(k/e)h²)であり、generic eのrank四は直接Pfから成立する。

## 5. 全pullback閉包への帰結

(1)のc'非零、Ω非零、Ω=c'零の全rank Hを尽くしたので、
**非零constant vについてi_vγ=0となる全最高部の全unimodular低次延長**は有限既知閉和へ入る。
γ零も含めると同じ帰結である。
condition存在vはproper projective incidenceで閉じるため、
S23全最高閉包、singular quadratic-pencilの全common-kernel、E3 rank一境界にも全点で適用出来る。
新成分零。閉log integrating-factor源の残例外、全h5分類の完成、scheme構造はここでは主張しない。
