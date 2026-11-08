# Ga9 linear-w は全 T36 の外：order二零点の二次jetと quadric-pencil の閉性

2026-11-08、PRIMARY E。純 whole-source／global-component 証明。
新しい point、jet/operator の計算、rank/minor、profile、GB 又は CAS は実行しない。
ここでいう二次jetは**既存全tensorの印刷式から手で取る固定の幾何不変量**であり、
局所変形jetの延長又は新anchor計算ではない。
結論は Ga9 unit dense23 orbit が whole closed T36 component に含まれないこと。
Ga9 の新成分性、次元、T37帰属／非帰属、全E完了は主張しない。

## 1. 四変数 quadratic bivector の proper quadric-pencil 条件

V を四次元vector spaceとする。quadratic homogeneous bivector Q に対し

```
alpha_Q = contraction(E_V wedge Q, volume_V)
```

は Euler に horizontal な cubic one-form、すなわち H0(P(V),Omega1(4)) の元となる。
volume 選択は非零 scalar だけを変え、GL(V) の作用と可換する。

Grassmannian Gr(2,Sym²V*) から

```
span(F2,G2) |-> [G2*dF2-F2*dG2]                            (1)
```

への map は全域で定義される。
実際、独立 F2,G2 で右辺が零なら、char零の有理関数体で d(F2/G2)=0 なので
F2/G2 は定数、独立性に反する。
基底変更で式は determinant 倍されるので Gr 上の regular projective map である。
Gr は proper だから、その image R22 は projective closed。
共通因子を持つ quadric pencil も除外せず、対応する非primitive cubic one-formをそのまま含む。
R22 の affine cone に零を加えた C22 は closed GL4-invariant cone である。

## 2. T36の全period source はこの jet 条件を持つ

source順は (A,B,C,u,v)、p=p0非零。
true completed source の vertex q=[1:0:0:0:0] は order二の projective zero である。
A=1 chartでは base W=(B,C,u,v) の全 quadratic part が

```
QBC=pBC+10t2t3uv/p,
QBu=3pBu/5+t2v², QBv=pBv/5,
QCu=-pCu/5, QCv=-3pCv/5+t3u²,
Quv=-puv/5.                                               (2)
```

vertex rows PAB=t0u²,PAC=t1v² による affine correction は cubic 以上なので、
二次jet(2)を変えない。これにより全t0/t1を捨てずに扱っている。

```
q+=(p/5)Bu+t2v², q-=(p/5)Cv-t3u²
```

と置く。全四 cubic coefficient の手展開で

```
alpha_Q = (5/p)*(q-*dq+ - q+*dq-).                        (3)
```

例えば dB係数は CQuv-uQCv+vQCu=u*q-、
dC係数は -(BQuv-uQBv+vQBu)=-v*q+、
du係数は BQCv-CQBv+vQBC=pBCv/5+t3Bu²+10t2t3uv²/p。
残りdv係数は -pBCu/5+t2Cv²-10t2t3u²v/p で、(3)の全式に一致する。
p非零なのでq+,q-は独立quadric。従って alpha_Q は C22 内。
period/edgeを固定した一つのsource点ではなく、全sourceの式(2)で成立する。

## 3. Global component への閉拡張

projective Poisson parameter space と P4 の積で、点qにおける sectionの0/1jetが零というincidenceを取る。
この条件下では quadratic leading tensor は座標変更に対し tangent GL4 と section scalar だけで変わる。
従ってその alpha_Q が C22 に入る条件は associated vector bundle 上の closed 条件。
alpha_Q=0 も含める。

q を忘れる射は P4 の properness により closed。
よって「ある order二以上のzero qがあり、alpha_Q が C22」という条件は
**全 projective tensor space 上で closed**。
T36 の一般 true sourceと全PGL像は§2でこの条件を満たすから、その image closureである
whole closed component T36 の全点も同じ条件を満たす。

normal-rootの二本のalignment、plane normal multiplicity、flat normalization、
constant preconjugation又は固定1PSの仮定を使わない。
従って moving source、singular basis、periodの任意のvaluationやzero-edge境界もこの閉性が覆う。

## 4. Ga9の全候補qはfiber lineで、全quadratic jetは非零

canonical四倍Ga9 tensorは

```
Pxy=4h²-y², Pxz=-hx+2yz+4w²,
Pxw=16hy-4wy, Pyz=hy, Pyw=0, Pzw=-4hw+4y²,
h row=0.
```

h非零なら Pyz=hy と Pxy からzero点は存在しない。
h零で E wedge P の yzw成分は4y³なので y=0、続くxzw成分は4w³なのでw=0。
逆にh=y=w=0では全tensorが零。
従って reduced projective zero support は厳密に

```
ell={h=y=w=0}=P(span(x,z)).
```

このlineの各qでは0/1jetが零になる。
base B0=(h,y,w)、fiber=(x,z)で mixed部は全てbase quadratic、fiber-pairのlinear-base項は
affine fiber chartではbaseとの積を一つ追加して二次以上となるためである。

固定D theoremの canonical Saito frame を

```
D=-4hF, F=4h³w-8h²y²+y⁴,
A=U+(4y/5)E_base, B=V+(4h/5)E_base,
det(E_base,A,B)=D
```

とする。Fはprimitiveなw一次式で絶対既約。
A,Bはquadratic logarithmic fieldsで、A(F)=-(4y/5)F、B(F)=-(4h/5)F。
q=[x:z]=[xi:zeta] における二次jetは、残るfiber tangent rとbaseを用いると

```
Qq=partial_r wedge Wq,
Wq = xi*B-zeta*A modulo (linear form)*E_base.              (4)
```

(4)は cotangent 全tensorをaffine fiber chartへ移して最低二次項を取るだけで従う。
fiber bracketの項は radial correction となり、alphaには寄与しない。
従って

```
alpha_Qq = contraction(E_base wedge (xi*B-zeta*A), volume_base)
```

を四変数base+rへpullbackしたものに非零scalarを掛けた形。
E_base,A,B はgeneric独立なので、どの [xi:zeta] に対しても alpha_Qq は非零。
しかも全qで F=0 の irreducible quartic coneがこのfoliationの invariant hypersurfaceとなる。
全ellを扱い、x非零chartの一点だけへ正規化していない。

## 5. Invariant irreducible quartic と quadric pencil は両立しない

非零 cubic one-form alpha が独立quadric pencil (1) から来ると仮定する。
その rational first integral は J=F2/G2。
一方§4のalphaは irreducible quartic hypersurface F4=0 を保つ。
alphaの全係数次数は三なので、F4が全係数を割ることはなく、F4一般点でfoliationは非零。
三次元P3のcodimension一foliationはそこで F4 の tangent planeそのものとなる。
従って dJ は F4 の全tangent方向で零、char零より J|F4 は定数c。
分母quadricはdegree四のF4を含めないので、その一般点で定義される。
ゆえに F2-cG2 が F4 上で恒等零。
左辺degree二、F4 irreducible degree四だから左辺は零polynomial、独立性に反する。

したがって alpha_Qq はどのq∈ellでも C22 に入らない。
Ga9は他にzero点を持たないので、§3のwhole T36 closed conditionを満たさず、
**Ga9 unit dense23 orbitは whole T36 の外**。

## 6. 射程

固定D normalizerでGa9 unit sourceを一actualorbitへ戻す既定理により、同じ全unit軌道に結論が及ぶ。
磁気zero軌道は別であり、そのT36非所属へこの結論を外挿しない。
T37 又は未知親のmembershipは残る。
T36 negativeから新component数、scheme smoothness、局所支持又は全E完了を推論しない。
本証明は全源閉性と既印刷式の手恒等式であり、新numerical certificateを仮装しない。
