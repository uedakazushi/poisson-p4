# 既知T37の零secantを含む全coisotropic plane：line+conicと三重因子の制限

2026-11-08、PRIMARY E、純手証明。
K4aのactual24親は全150/PGL/Euler恒等式と独立reviewにより既知T37へ認識済み。
以下はその同じ整数generic tensorの全secant-plane classを扱う。
stable二源の全帰属又は全T37非所属を主張しない。

## 全coisotropic plane曲線

`K4A_SMOOTH_PARENT_ZERO_SECANT.md` の actual zero line
`ell: a=3c,u=b,v=9c` を使い、normalを

```
z=a-3c, w=u-b, h=v-9c
```

とする。actual全二次tensorからの三 normal bracketは正確に

```
{z,w}=810cz+60bw+258z²+18w²,
{z,h}=60bh-72zw+36wh,
{w,h}=-810ch-288z²-18w²-54zh-36h².                          (1)
```

Euler ambiguity `E∧L`、`Lz=-60b,Lw=810c,Lh=0` を除けば、
linear-along-ell項は全て消え、残るnormal bracketは(1)のquadratic項だけである。
これはactual Euler representative変更であり、plane coisotropyを変えない。

ellを含む任意のplaneはprojective normal direction n=(z0,w0,h0)∈P2で表される。
coisotropyは n と normal bivector dualとのpairingであり、全条件は一つのcubic

```
-6*(h0-3z0)*(3w0²+6z0*h0-16z0²)=0.                        (2)
```

したがって全plane familyはlineとsmooth conicの和である。
conicはθ²=-6として

```
n(s,t)=(3s²,-θ*s*t,8s²+t²)
```

と全域でparameter化できる。accepted T37 frameのquarticをMöbius変更すると

```
Q(s,t)=[3s⁴+s²t² : 6s⁴+6s²t² : -3θs³t+θst³ : -2θs³t : 19s⁴+16s²t²+t⁴].
```

そのnormal三座標は厳密に `(t²-s²)*n(s,t)` なので、これはQ(s,t)とellをspanしたplane family。
端点t=±sでは同じ式のproperなplane極限を取る。
lineとの交点は t=±s、すなわち `(z0,w0,h0)=(3,±θ,9)` の二点。
このFano schemeの特別flagをgeneric smoothflagとして扱わない。

## Actual quintic は ellに沿ってorder三以上

affine c=1、line coordinate b と上記normalを取る。
零line上のactual Poisson firstjetは

```
pi1=∂b ∧ X1(z,w,h;b)
```

であり、normal-normal firstjetは全て零。
これは同noteの六affine entryの全firstjet（pi_au=pi_ab、pi_av=0、pi_uv=pi_bv）から直接従う。
従ってpi1∧pi1=0、actual Pfaffianのnormal二次項は零。
元tensorもell上で零なので、Dはellに沿ってorder≥3である。
線形idealではsymbolic/ordinary powersが一致し、normal係数はline上のpolynomialだから

```
D ∈ I_ell³
```

が全線上で成り立つ。別点や追加operatorを使わない。

どのplane Sがellを含んでも、sectionは `D|S=L_ell³*Q2` 又は零。
(2) のline/conic両族と二intersection flagsを全て含む。
先のspecial coisoplaneの `c²(v-9c)³` はこの一般条件の一例である。

## Uniform normalcone classの正確な停止条件

coisotropic planeに二normal方向を一様にepsilon倍し、tensorをepsilon倍すれば、
rank四の正則normalconeのactual quinticは元のplane sectionそのもの。
regularなplane移動やそのmagnetic defectを許しても、正則base frameを保つこのclassの端点divisorは
同じsectionであり、三重因子は失われない。
cotangent表示の全magnetic coefficientはanchor determinantを変えないこととも一致する。

従ってこの **ellを含むplaneのuniform normalcone class** では、
二重Hと他のsimple因子を持つ stable_line_conic/stable_three_lines も、
二つのsimpleHと既約cubicの平面残三源も得られない。
Fano nodeからessentialmagnetic方向が生じる可能性だけでこのdivisor障害は解消しない。

これは既知T37 component全体の陰性ではない。
ellを含まないcoisotropicplanes、singularbaseframe、別のweightsやactualmoving-source familyは残る。
同じuniform-plane profileを増やさず、次のwhole-source幾何へ切替える。
