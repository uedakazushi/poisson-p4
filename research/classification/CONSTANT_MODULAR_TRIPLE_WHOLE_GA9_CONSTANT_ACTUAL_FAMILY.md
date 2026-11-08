# constant TRIPLE の全源は Ga9 非零定数 c の既23成分へ入る：実Poisson多項式曲線

2026-11-08、PRIMARY E。担当 `/root/plane_collision_review`。
**新純証明producer、独立査読待ち。** 新しいoperator、point、rank/minor、jet又はGB無し。
同じforcedcanonicalからの一つのactual ambient曲線を全六括弧で明示する。
非線形Darboux写像はPoisson曲線の構成・全Jacobi証明にだけ使い、PGL所属の代わりには使わない。
一般fiberの既知成分帰属は、別の**実際の可逆affine/PGL frame**と全source定理で証明する。

## 1. 受理済み全canonical source

物理座標(t,x,y,r)でforced TRIPLE representativeは

```
Ptx=x/3+y, Pty=−xy+y/6−r/3,
Ptr=−3y²+xr−5y/4+r/2,
Pxy=1, Pxr=−3x+3/2, Pyr=3x²+3y.
```

Pfは f0=r+H、H=x³+3xy−3y/2、divP=−∂t。
全最高triple-line sourceは、既受理の必要十分finite reductionによりknown閉和又は
この一つのactual orbitに入り、残余全四parameter源は分母を追加しないactualnormalizerで
**全域**同一orbitへ移る。
本稿がこのrepresentativeを既知閉成分へ送れば全残余source、従って全highest-triple枝が閉じる。
この全source→orbit橋は別入力で、単なる一標本のJacobi成功から推測しない。

## 2. Poisson productと全parameter可逆な構成

補助空間の座標を(T,x,y,s)とし、元のPoisson bracketを

```
{x,y}=1, {T,s}=s, other pairs=0
```

と取る。二つの二次元Poisson構造の直積なので全Jacobiが零、signed Pfはs。
以下を置く。

```
a=x/3,
k=xy/3+y²/2−x⁴/12,
A=e/12, B=e*(15−e)/72,
M=A*x²*y+B*y².
```

parameter eについて次のpolynomial automorphismを取る。

```
r=s+e*T−H,
t=(1−e*a)*T−a*s−k+M,
x'=x, y'=y.                                             (1)
```

(x,y)固定の(T,s) blockのJacobian determinantは
(1−ea)*1−(−a)*e=1。全四変数Jacobianも一。
全e、全点で逆は正確に

```
T=t+a*(r+H)+k−M,
s=r+H−e*T.                                              (2)
```

であり、1−eaを反転しない。
従って(1)はeの全域で可逆なpolynomial写像。
e=0の補助Darboux逆は t=T−k−a*s、r=s−Hで、§1の全六括弧へ戻る。

## 3. 全六physical括弧はquadratic polynomialである

(1)で押し出したbracketをP_eとする。全四変数・全eでPoissonであり、
直接chain ruleから

```
Pxy=1,
Pxr=−3x+3/2,
Pyr=3x²+3y,
Ptx=x/3+(1−2B)*y−A*x²,
Pty=(−1+2A)*x*y+y/6−r/3,
Ptr=−e*t+(1−e/3)*x*r
    +(-3+3*e/4+e²/8−e³/72)*y²
    +(5*e/12)*x*y−5*y/4+r/2.                            (3)
```

高い次数を捨てていないことを明示する。
{M,H}=M_x H_y−M_y H_x であり、

```
{M,H}=−3A*x⁴+(3A−6B)*x²*y−3A*x*y−6B*y².
```

K=k+a(r+H)=x⁴/4+x²y−xy/6+y²/2+xr/3 と書くと

```
Ptr=(Ptr)_0−e*t−e*K+e*M−{M,H}.
```

x⁴係数は−e/4+3A=0、x²y係数は−e+eA−3A+6B=0。
残りのxy,y²,xr係数が正確に(3)になる。
他五括弧にも三次以上は無い。
従って(3)は**全parameterで元の二次Poisson係数空間に入る実多項式族**である。
元productの全Jacobiと可逆pushforwardから、全physical四Jacobi及び斉次化後の元350projective Jacobiは
全eで恒等零。有限jetの延長を仮定していない。

## 4. 実Pfaffian、modular場、全zero補集合

det(1)=1だから、全signed Pfは(2)のsそのものである。

```
f_e=r+H−e*[t+a*(r+H)+k−M].                              (4)
```

これは(3)からのPf=Ptx*Pyr−Pty*Pxr+Ptr*Pxyとも一致する。
e=0でf0、e≠0でx⁴係数は−e/4、t係数は−eである。
よって全e≠0でf_eはtについてleading constant unitの一次式、従って絶対既約。
四次斉次化F_eもh因子を持たず、絶対既約quartic。
actual projective divisorは

```
e≠0: D_e=h*F_e (simple H plus irreducible quartic),
e=0: D_0=h²*(x³+3xyh−3yh²/2+rh²).
```

点(t,x,y,r)=(0,0,0,1)では全eでf_e=1なので、族のprojective係数又はrank四を失うparameterは無い。

volume Jacobian一によるmodular共変性、又は(3)の全偏微分から

```
Z_e=(-1+e*x/3)*∂t−e*∂r.                                (5)
```

特にe≠0の線形部はrank一nilpotent、定数部はそのimageの外でshort型。
**ここでPfQ0を仮定していない。** 実際e≠0ではdegree四Pfが非零なので、
multiple-Hの全17源定理を使うと仮定違反になる。
次節はsimple-H＋既約quarticの正しいGa9入口を使う。

## 5. generic全fiberのactual affine Ga9 frame

e≠0に対して元座標を、次の新座標(X,Y,Z,W)へ送る。

```
X=t−r/e, Y=e*x/3, Z=−r/e, W=y.                         (6)
```

これは実際のaffine GL4、det=1/3、tensor scale一、Euler補正零。
逆も全て一次式で

```
r=−e*Z, x=3Y/e, y=W, t=X−Z.
```

(5)は正確に

```
Z_e=Y*∂X+∂Z
```

へ移る。これは全Ga9短translation cofactorの固定標準形そのもの。
Ga9 invariant u=X−Y Z は

```
u=t−r/e+x*r/3.
```

(4)は

```
f_e=−e*u+G_e(x,y),
G_e=(1−e*x/3)*H−e*k+e*M.                               (7)
```

と書ける。
(6)のvolume determinant一/三を戻した新signed Pfはf_e/3。
Ga9全sourceのstraightening s_flow=Z、u=X−YZ におけるpotential Φ は
Pf=−Φ+constant の規約だから、

```
Φ_u=e/3 ≠0 (constant).                                 (8)
```

従って全e≠0 fiberは、**全係数を保つGa9 nonzero-constant-c source**に入る。
これは元28parameter全源のnecessary/sufficient分類からの帰結であり、
potentialの形が似ていることだけによる親推測ではない。
(3)の二次次数、全Jacobi、divergenceの固定、既約quarticは全て上で確認した。
Ga9 c=1へ規格化する非零constant尺度も元の全source定理のactual affine/scalar操作である。
quartic leading coefficientが全e≠0で非零なので、そのfullsourceのr非零というdegree四条件も失わない。

既独立受理済み `GA9_CONSTANT_COMPONENT.md` と `GA9_CONSTANT_ANCHOR_REVIEW.md` により、
この完全sourceの全PGL像は既23成分

```
R4_GA9_CONSTANT_HYPERPLANE_QUARTIC_23
```

に入る。一つのanchorだけの局所結論を全sourceへ拡張したのではなく、
元定理で証明済みの全十源包含を使う。

## 6. whole TRIPLEと全constant-modular枝への結論

(3)の全係数はeのpolynomialであり、e=0は指定されたforcedcanonicalそのもの。
e≠0は同一の既知closed Ga9 c1 componentに入るので、閉性によりe=0も入る。
このcomponentはPGL/scalar不変だから、全canonical orbitとその閉包も含む。
§1の全source inverseを戻し、**最高triple-lineの全constant-modular源は既知閉和K**へ入る。
新成分増分零。元のwhole4sourceの個々のp,d,l,Lやunit補集合を落としていない。

これは新しいactual正包含であり、直前のQ51/Q52全graph陰性と矛盾しない。
旧TRIPLE tangent28/orbit22というnonclosure reportもそのまま正しい。
DOUBLEの親認識又は新成分性、全E完成はこの証明からは決めない。
