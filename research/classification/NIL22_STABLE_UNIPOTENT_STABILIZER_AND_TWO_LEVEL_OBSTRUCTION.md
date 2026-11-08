# Stable 二源の全 connected stabilizer と二段階 gap一の実変形障害

2026-11-08、PRIMARY E。純手証明、新CAS・operator・anchor・profileを使わない。
whole-source の既全磁気商と actual divisor から、uniform normalcone 反復の構造的限界を示す。
全 MPS 親の非所属又は stable の成分性は主張しない。

## 1. intrinsic base quotient の connected stabilizer は自明

stable 二源の actual divisor は、三座標だけに本質的に依存する

```
line_conic: h²*L*Q, Q smooth, L tangent, h secant, h∩L off Q;
three_lines: h²*L1*L2*L3, 四線に三重交点無し。
```

全五 partial の共通 annihilator が二次元 vertex V を intrinsic に定める。
従って projective tensor stabilizer は V と商 P2 を保つ。
line_conic では Q、二重 h、単純 L が各々保存される。
Q の normalization P1 上、h の二交点と L の接点は相異なる三点である。
これらを保存する PGL2 は有限なので、base P2 の connected stabilizer は自明。
three_lines でも一般位置の四直線の保存群は有限なので同じ結論。

base を恒等写像に正規化した connected projective stabilizer の全候補は

```
base'=base, q'=A*q+ell(base), A∈GL2, ell∈Hom(C3,C2).
```

これは vertex 保存からの全 block 行列であり、特定 normalizer 部分群だけではない。

## 2. 全 sixmag quotient による scalar の除外

元の全 source は h=1 で

```
{qi,pj}=Mij(p), {p1,p2}=0,
{q1,q2}=R(p)+linear_in_q,
R=m0+m1*x+m2*y+m3*x²+m4*xy+m5*y².
```

二 anchor は C 上線形独立である。
q'=Aq+ell の mixed brackets が projective scalar λ 倍の元 tensor と等しければ
`A*M=λ*M`、従って A=λI2。
q'=λq+ell の全 q,q bracket を projective scalar λ で割ると、磁気 class は

```
[R] -> λ*[R] in sixmag / Im(delta).
```

ell の作用は既受理の全六変数 affine translation delta そのものである。
line_conic で rank delta=4、three_lines で rank delta=3。
その全逆と二・三残不変量 I は `NIL22_STABLE_PROJECTIVE_ACTION.md` に保存済み。

I≠0、すなわち residual magnetic class が非零なら、固定条件から λ=1。
残る全 stabilizer は厳密に `ker delta` の translations である。
したがって connected projective stabilizer はそれぞれ

```
stable_line_conic: Ga²;
stable_three_lines: Ga³.
```

両群は unipotent で、非自明な代数的 C* subgroup を持たない。
これは全 I≠0 開部分について成立し、特定の新 witness に依存しない。
I=0 の補集合は残す。このとき fiber scaling C* は実際に存在し、
three_lines の同境界は既全 T40 graph 認識とも整合する。

## 3. gap一の全 regular moving frame は pure limit と同じ orbit

任意の finite-dimensional PGL 表現を取り、二段階の cocharacter を

```
G_e=diag(I_r,e^(-1)*I_s), r+s=5, r,s>0
```

とする。regular invertible frame B(e) を B(0) で分け、まず B(0)=I としてよい。
全 block を展開すると

```
G_e B(e) G_e^-1
 = [[I+O(e), O(e²)], [B21^(1)+O(e), I+O(e)]].
```

右辺は e=0 でも正則可逆、極限は lower unipotent。
従って `e^m*(G_e B(e))_*P` の非零極限があるなら、それは
`e^m*(G_e)_*P` の pure cocharacter limit の正則 linear conjugate である。
先頭に任意の regular invertible A(e) を掛けても、その A(0) conjugacy が加わるだけ。
pure limit の projective class は非自明な G(C*) に固定されるので、
全 regular moving-frame limit も conjugate C* subgroup に固定される。
この線形群の等式は formal arc 全次数に成立し、有限 jet の外挿ではない。

## 4. 全 uniform two-normal contraction の正確な限界

特に r=3,s=2 の一様 normalcone と、全 regular pre/post frame を含む class は
I≠0 の stable 二源のどの点にも到達しない。
first-order plane-motion defect が anchor translation coboundary に留まることは
この全行列等式の微分版であり、欠けた residual magnetic class を自由に加えられない。

これは全 coisotropic plane を仮定せず、仮に全 plane が分類されていなくても使える。
section が目的配置に一致しても、この stabilizer 障害は消えない。

全親の非所属へは上げない。二つ以上の weight gap、singular moving preframe、
又は orbit 外の genuine source-parameter deformation は本補題の範囲外。
別の weight を無根拠に走査することも正当化しない。
次の正の橋は full residual magnetic class を変え得る数学的構造を先に備える必要がある。
