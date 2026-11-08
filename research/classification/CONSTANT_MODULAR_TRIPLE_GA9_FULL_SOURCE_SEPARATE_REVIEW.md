# constant TRIPLEの全Ga9定数源への実曲線：独立純査読

2026-11-08、PRIMARY E_multiple_hyperplane、担当 `/root/exceptional_coverage`。
**独立手展開・旧全源定理の仮定照合。新算術実行無し。結論PASS、数学的欠陥無し。**
rootによる別exact universal再構成0.315秒は、その証明書readbackとして区別する。

## 1. 全parameter polynomial inverseと六physical括弧

H=x³+3xy−3y/2、a=x/3、k=xy/3+y²/2−x⁴/12、
A=e/12、B=e(15−e)/72、M=A x²y+B y²とする。
補助productは{x,y}=1,{T,s}=s、他零でPoisson。

```
r=s+eT−H,
t=(1−ea)T−as−k+M.
```

逆T=t+a(r+H)+k−M,s=r+H−eTを両向きに代入すると恒等となる。
(T,s)block determinantは(1−ea)+ae=1、全Jacobianも一。
1−ex/3を反転せず、e=0を含む全域でpolynomial automorphismである。
この非線形補助写像をPGL所属の代用にはしない。

手でchain ruleを展開した。Pxy=1、Pxr=−3x+3/2、Pyr=3x²+3y、
Ptx=x/3+(1−2B)y−Ax²、Pty=(−1+2A)xy+y/6−r/3。
最後のPtrではK=k+a(r+H)=x⁴/4+x²y−xy/6+y²/2+xr/3を用い、

```
Ptr=(Ptr)_0−et−eK+eM−{M,H},
{M,H}=−3Ax⁴+(3A−6B)x²y−3Axy−6By².
```

x⁴係数−e/4+3A、x²y係数−e+eA−3A+6Bは厳密に零。
残る全式は

```
Ptr=−et+(1−e/3)xr+C_e*y²+(5e/12)xy−5y/4+r/2,
C_e=−3+3e/4+e²/8−e³/72.
```

したがって全eの全六括弧が二次以下である。
productのJacobiと可逆pushから、これは有限jetでなく実際のPoisson polynomial曲線。
e=0は指定されたTRIPLE全六括弧と一致する。
rootの別exact reportは全150 canonical、全350、全126 actualDをこの同じ曲線でPASSしている。
その全行列を本査読で再実行したとはしない。

## 2. 実Pfaffian・cofactor・全zero補集合

Jacobian一だからsigned Pfはs=r+H−e[t+a(r+H)+k−M]。
e≠0ではt係数−eが定数単元で、x⁴係数−e/4も非零。
従って絶対既約quarticで、斉次化はhを因子に持たない。
実Dは単純Hと既約quarticの積。e=0だけでh²Cへ退化する。
全eで(t,x,y,r)=(0,0,0,1)のPfは一なので、曲線のtensorは零にならない。
発散の全四係数を手で微分すると

```
Z_e=(−1+ex/3)∂t−e∂r.
```

これはsimple-HのGa9入口であり、Pf(Q)=0を仮定するmultiple-H十七源定理を適用しない。

## 3. 実際のlinear frameとc=1への全逆

全e≠0について

```
X=t−r/e, Y=ex/3, Z=−r/e, w=3y/e.
```

これはdet=1/eの実GL4変換、tensor尺度一、Euler補正零である。
逆はt=X−Z,x=3Y/e,r=−eZ,y=ew/3。
Z_e(X)=Y、Z_e(Z)=1、Z_e(Y)=Z_e(w)=0なのでmodularは厳密にY∂X+∂Z。
最後のw尺度は中間frame W=yにおけるc=e/3をw=W/cにする操作そのもので、
{Y,w}=1となる。signed Pfも中間Pfをcで割り、最終Pf=s/e。
u=X−YZの係数が−1だからPhi_u=1である。

この実GL4 pushを全六括弧で手展開すると、b=3/e、r_s=81/e⁴として

```
P_YZ=bY−1/2, P_Yw=1, P_Zw=r_s Y²+bw,
P_XY=1/2+(1/3−b)Y−Y²/4+2F2*w,
P_Xw=Z+r_sY²+(1/2−b)Yw+(b+1/6)w,
P_XZ=X−Z/2+(b−1)YZ+A2,
A2=F2(1+2b)w²−5Yw/12+5w/12,
F2=e²(e−3)(e−12)/648.
```

最後のw²係数では−e*C_e/9=F2(1+2b)を展開確認した。
したがってpotentialが似ているだけの認識ではなく、全六physical係数のliteral一致である。

## 4. 既necessary/sufficient全十源へliteralに入る

GA9_CONSTANT_COMPONENT §1、元HYPERPLANE_QUARTIC_NIL211_SHORT_TRANSLATION §7と
その独立査読§4を読み、全28源からの必要十分性と全仮定を照合した。
該当十parameterは

```
r=r_s=81/e⁴, b=3/e, alpha=−1/2,
p=q=D=G=k=0, A=−27/e⁴=−r_s/3, H=1/2.
```

元定義B_s=b−1/4、E_s=alpha−[(1−3b)A+qB_s]/r_sからE_s=−b−1/6。
元定義F2=−(1−b)B_s/(2r_s)も上のe式と一致する。
Vy=1/2−bY、Vw=bw+r_sY²、

```
phi=r_sY⁴/4+A Y³+B_sY²w+E_sYw+F2w²+H w,
Phi=u+phi, Vu=−u−A2.
```

A2=phi+Vy*phiY+Vw*phiwを独立に項別展開した。
Y⁴係数r_s(1/4−b+B_s)、Y³係数A(1−3b)+r_s(E_s+1/2)、
Y²w係数B_s(1−b)+2r_sF2、Y²係数3A/2+r_sHは全て零。
残る係数はE_s+B_s=−5/12、H+E_s/2+bH=5/12、
w²係数F2(1+2b)のみで、§3のA2そのものになる。
元変数への逆u=X−YZ,Z=flowの全六括弧も§3に一致する。

r_sは全e≠0で非零であり、全十源が要求する唯一のquartic openを満たす。
e=3又は12でF2=0となるが、元定理はF2又はb−1等を反転していないので両値を含む。
e≠0の他の値・低次parameter零も全て同じnecessary/sufficient源に含まれる。
全二次係数、全Jacobi、固定modular、単純H＋既約quarticを確認済みなので、
未表示の線形源又はpotential compatibilityを仮定していない。

## 5. 全四sourceとcomponent閉性への橋

既独立受理されたGa9非零定数c全十源の飽和像閉包は、既23成分
R4_GA9_CONSTANT_HYPERPLANE_QUARTIC_23そのもの。
従って全e≠0のP_eは同じ既closed componentに入る。
P_eの元二次係数はeのpolynomialで、e=0の係数も零tensorでないから、閉性によりP_0も入る。

受理済みTRIPLE全四sourceは、全parameterで実逆を持つactual normalizerにより一軌道へ移る。
既成分はPGL/scalar安定なので、その全軌道と閉包、従って全四parameter源も含まれる。
元全最高triple-line還元の他のknown枝と合わせて、この最後の残余を閉じる。
一点のsmoothnessを境界から逆推論したり、potentialの類似だけで親を推測したりしていない。
新成分増分は零。旧T39/T40 graph陰性と両立し、そのbytesは変更しない。

**このactualfamilyと全Ga9十源への実所属を独立純査読PASSとして受理する。**
元producer、root算術、共有code/job、台帳又はGitは本査読で変更していない。
