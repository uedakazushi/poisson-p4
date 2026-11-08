# constant TRIPLE：whole T01/T21/T36/T38/T40 と Ga9 の純幾何排除

2026-11-08。PRIMARY E。担当 `/root/plane_collision_review`。
**producer pure proof、独立査読待ち。** 新point/operator/minor/局所変形jet/GB無し。
固定 actual tensor の全零台と最低二次項を手で展開し、proper closed 全成分条件と比べる。
親の正同定、新成分性、全E完成は主張しない。

## 1. 全actual representativeと因子

physical座標(h,t,x,y,r)、h-row零、六括弧は

```
Ptx=xh/3+yh,
Pty=−xy+yh/6−rh/3,
Ptr=−3y²+xr−5yh/4+rh/2,
Pxy=h², Pxr=−3xh+3h²/2, Pyr=3x²+3yh.
```

これは全highest-triple必要十分源の一actual orbit正常形である。
integer ambient witness は固定diag(12t,x,2y,2r/3)の同orbitであって、別点ではない。
実退化因子は

```
D=h²*C, C=x³+3xyh−3yh²/2+rh².
```

Cはrについてprimitive linearで、leading h²と定数部のgcdは一。
よって絶対既約であり、unique double-Hはh。
Cのessential rankは四である：constant方向D0がCをannihilateするなら、
rh項からD0h=0、yh項からD0x=0、xh項からD0y=0、最後にD0r=0。
Dそのものも四変数を本質的に使う。

## 2. 全intrinsic zero台と全order二候補

h=1ではPxy=1なのでprojective零点は無い。
h=0における全Euler wedge四cubic、順序(txy,txr,tyr,xyr)は

```
x²y, x(3y²−xr), 3tx²−2xyr+3y³, 3x³.
```

最後からx=0、続いて第三からy=0。
逆にh=x=y=0では全十Euler cubicsが零である。
従って厳密に

```
Zero(P)_red=L={h=x=y=0}≅P1_(t:r).
```

r=1 chartのaffine bracketはPij+xi*Pjr−xj*Pir。
L上の任意の有限点t=t0でπtyのlinear termは−h/3。
従って全r≠0点ではsectionのfirst jetは非零、order二零点ではない。
唯一残るq∞=[h:t:x:y:r]=[0:1:0:0:0]ではt=1 chartを使う。
全0/1jetは零で、二次最低項は

```
Q=h²∂x∧∂y+(−3xh+3h²/2)∂x∧∂r+(3x²+3yh)∂y∧∂r.
```

h-rowは二次最低項零であり、これは全六affine chart entryを残して得る。
全order二零候補はこの一点だけである。

## 3. T36のproper quadric-pencil条件と合わない

t=1 tangent空間のEulerをE4、volumeをdh∧dx∧dy∧drとすると、
最低項から得る cubic Euler-horizontal one-formは符号を固定して

```
α=i_(E4 wedge Q)vol4=3C*dh−h*dC.
```

四coefficientは

```
dh: 3x³+6xyh−3yh²/2+rh²,
dx: −3hx²−3yh²,
dy: −3xh²+3h³/2,
dr: −h³.
```

特にαは非零。αに対応するfoliationはirreducible cubic hypersurface C=0を保ち、
その一般点で非零である。
後者は、仮にCが全 cubic coefficientを割ればα=C*constant one-formとなり、
Euler horizontalityから零を強制されることからも分かる。

全truecompletedT36は「あるorder二zeroのαがquadric pencil <F2,G2> の
G2dF2−F2dG2のclosed affine cone（零含む）に入る」というproper closed条件を持つ。
全source四coefficient恒等式、Grassのproperness、全marked点を忘れるpropernessは
既別手照合 `CONSTANT_MODULAR_DOUBLE_WHOLE_T36_T38_EXCLUSION_SEPARATE_REVIEW.md` §4を使う。
したがってmoving frame/edge/periodの任意境界も同じ条件を満たす。

TRIPLEの唯一候補q∞でこの条件が成立すれば、F2/G2がrational first integral。
C一般点でdefinedかつ全tangent方向で微分零なので、C上で定数c。
よってdegree二のF2−cG2がirreducible degree三Cを割り、独立性に反する。
従って **TRIPLEはwhole T36の外**。

## 4. T38/T40のwhole zero-surface条件

truecompletedT38全generic sourceには、既proofのexplicit twisted-cubic cone zero surfaceがある。
T40=q52全generic sourceには、既受理のelliptic rank二degree五bundle W上の
smooth elliptic scroll S0=P(W)がrank零surfaceとして存在する。
これは `RANK_ZERO_QUINTIC_SURFACE_GENUS_ONE_RECOGNITION.md` §6 の
全q52 normalization/conductor入力と既受理MPS/q52全構成を使う。
一つのunstable特殊sourceがT40に入ることから逆向きにgeneric geometryを推測していない。

各全generic family/PGL imageはそれぞれ全T38/T40に稠密である。
proper universal projective zero incidenceのfiber dimension上半連続性より、
「dimZero≥2」は両whole閉成分に延長する。
TRIPLEの全zero台は§2で一次元なので、**whole T38/T40の外**。

## 5. T01/T21/Ga9のwhole factor/essential-variable条件

T01の全generic sourceのDは五つのlinear factors。
五fold projective multiplicationのproper imageはclosedなのでwhole T01も全D all-linear。
TRIPLEはirreducible cubic因子を持つからその外。

T21の正しい全sourceでは

```
D=x3*x4*[p0*(t2*x0³+p1*x0*x1*x2−t1*x1³+t0*x2³)].
```

従ってDは「二linear factorsと、三変数以下のcubic」というproper closed incidenceに入る。
三変数のsubspaceをGr(3,V*)上で動かし、projective cubic/factorと共にproperに忘れれば、
whole T21の全境界も含まれる。
TRIPLEの唯一prime cubic C はessential四なので、この因子型へは入らない。
従ってwhole T21の外。
DOUBLEのcubicはessential三なので、DOUBLEへこの議論を流用せず、既torus orbit上界を使う。

Ga9 cuspidal quartic24成分の全generic Dは三変数で書ける。
全D essential rank≤3のlocus自体がGrass incidenceのproper imageでclosed。
TRIPLEのDはessential四だから、**whole Ga9成分の外**。
この排除はDOUBLEのD5 miniversal証明を必要としない。

## 6. 範囲

本稿が新たに排除したwhole parentsはT01,T21,T36,T38,T40,Ga9。
全long25の別properflag証明、他のfactor閉条件とは独立にscopeを保つ。
T39、T37、他のより小さい成分・未知支持の帰属/相違はこの稿からは決まらない。
projective tangent28を根拠に25親を排除したり、zeroのscheme厚みを計算したりしていない。
