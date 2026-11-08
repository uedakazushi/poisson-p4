# 共通factor付き三次theta：全低次fiberを保つ対数incidenceへの直接入口

2026-11-08、PRIMARY `E_multiple_hyperplane`。producer純証明、独立/root査読前。
`UNIMODULAR_H5_SINGLE_ACTUAL_ORBIT_FRONTIER.md` の最高部入口を独立に補完する。
旧G証明、旧六閉包記録、旧immutable reportは変更しない。新算術無し。

## 1. θfactorとγfactorを区別する

γはclosed homogeneous quadratic二形式、γ²=0、θ=i_Eγ/4、dθ=γ。
従ってi_Eθ=0、θ∧dθ=0である。
γの係数に共通factorが無くても、θに共通factorがあることは有り得る。
本稿は後者をγfactorの定理へ無断転用しない。

**補題。** 非零integrable cubic one-form θがEuler-horizontalで、その係数がnonconstant gcdを持つなら、
非零quartic D（two linear factors×one quadratic、repetition可）が存在して

```
D*dθ−dD∧θ=0.                                      (1)
```

を満たす。すなわちθは全closed logarithmic integrating-factor incidenceに入る。
この結論の後、各最高θに対する**全C6/L24 extension**の新pole定理を適用する。
最高θのgeneric densityから低次fiberを平坦と仮定しない。

## 2. gcd次数二

θ=F2 η、ηはhomogeneous linear一形式とする。
F2η∧d(F2η)=F2²η∧dηだからηはintegrable。
i_Eη=0よりη=i_E(dη)/2。
dηはconstant二形式で、その平方はd(η∧dη)=0。
η非零だからrank(dη)=2で、二つの独立linear coordinates x,yを取り

```
η=x dy−y dx
```

と出来る。D=x*y*F2は非零quarticで
θ/D=dlogy−dlogxがclosed。F2は任意quadratic、全rank・全factorを保持する。

## 3. gcd次数一：dηの二ruling

θ=F1 η、ηはhomogeneous quadratic一形式。
再びη∧dη=0、i_Eη=0、η=i_E(dη)/3。
ξ=dηはcoefficients linearのclosed二形式であり、
ξ²=d(η∧dη)=0。
従ってξ(V)⊂Λ²V*は線形完全等方空間である。
decomposable二形式線形空間の初等二ruling分類から次の全二場合へ分かれる。

### 3.1 共通constant一形式

ξ=dx∧b1、b1はlinear一形式と出来る。
dξ=0からdx∧db1=0であるため、全b1のx以外の三成分はquadratic G2の偏微分として積分できる。
従ってξ=dx∧dG2、G2はnonzero homogeneous quadratic。
Euler式から

```
η=(x dG2−2G2 dx)/3,
D=F1*x*G2,
θ/D=(dlogG2−2dlogx)/3.                              (2)
```

よって(1)。F1がxに比例する場合もD=x²G2となって同じ(2)が成立する。
G2が可約、singular、F1とfactorを共有する場合も除外しない。
G2がxのみの式ならξ=0、η=0となり元の非零仮定に反するので、余分なunit条件は要らない。

### 3.2 三変数support

ある三次元constant一形式空間Wに対しξ∈Λ²W。
Wを消すnonzero constant vについてi_vξ=0、closednessからLie_vξ=0。
η=i_Eξ/3もi_vη=0、Lie_vη=0となり、η自身が三変数pullbackである。

そのbaseを(x,y,z)とすると、Euler-horizontal quadratic一形式の全空間は八次元で、

```
A∈sl3 ↦ η_A=i_E i_(Ax)(dx∧dy∧dz)
```

が線形同型である。kernelはAxがxに比例するscalar行列だけであり、trace零では零。
targetの八次元はEuler写像の係数比較、又はΩ1_P2(3)のEuler列から得る。
この全linear parameter空間を使い、generic semisimple Aだけを全sourceと扱わない。

相異なる三固有値を持つAは非空稠密開で、actual base GL3のeigenframeにより

```
η=xyz*(λ dlogx+μ dlogy+ν dlogz), λ+μ+ν=0.
```

従ってD=F1*x*y*zに対しθ/Dはclosed（F1のresidueは零）。
F1はbase外の第四coordinateにも依存してよく、ηpullbackからγpullbackを推論していない。

quartic Dがfour linear factorsの積であるprojective parameter集合はproperであり、
(1)のθ/D incidenceをθへ射影した像は閉である。
固定F1を掛けた全η_Aの稠密semisimple開がこの閉像へ入るので、全Aも入る。
全Jordan、repeated eigenvalues、零residue、dependent F1もこの包含に残る。
ここで示すのは最高部の**閉incidenceへの包含だけ**であり、未知の低次extensionの極限存在ではない。
その後に使うpole全源定理は、incidenceの各個別θに対して全C/Lを残したものである。

## 4. 全非saturated入口と分類への接続

gcd次数三ならηはconstant、i_Eη=0からη=0なので非零θでは不可能。
従って§§2–3で全non-saturated cubic θを尽くした。
θ零はγ=dθ零であり、既全零最高部CAへ戻す。

θのgcdが一なら特異zero locusはcodimension二以上であり、通常のdegree二foliationの
六閉包定理をその本来のsaturated scopeで適用出来る。
共通factor付きθには本稿を適用して全quartic integrating-factor incidenceへ入れ、
`UNIMODULAR_LOG_POLE_WHOLE_BETA_ZERO_AND_SINGLE_ORBIT.md` の全extension定理により
有限既知閉和又は残一軌道P_*へ送る。

したがって新h5 aggregateの全入口は、非saturated locusについてLPT拡張の解釈へ依存せず閉じる。
外部六閉包定理の再証明、G非優越の再開、全h5残軌道の既知親同定はここでは行わない。
