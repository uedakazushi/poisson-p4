# constant modular：最高double-line+lineの全sourceを一軌道閉包へ有限還元

2026-11-08、PRIMARY `E_multiple_hyperplane`、担当 `/root/plane_collision_review`。
**純全係数証明。新operator・点・minor・kernel計算・profile・jet・GB無し。**
対象はA=0,c≠0、f3が相異なるdouble line+simple lineである全source。
全源を既知閉和と、最後に明示する**一つのactual orbitの閉包**へ還元する。
その最後の軌道の既知親／成分性は未解決であり、有限化だけでEを閉じたとは称さない。
producerと既依存pure証明の独立受理は別途必要である。

## 1. 全入口と二つの既知十分枝

κ=1、base順序(x,y,r)としてP=∂t∧Y+Jac(f)、divY=1、Yf=fを使う。
actual linear frameとf3に現れないrの非零尺度によりf3=x²yとする。
f3の全体係数を消すr尺度は可逆であり、f3=0は入口の外。
全Y2/Y1/Y0、f3/f2/f1/f0を最後まで残す。

以下の既知十分枝を使う。

* Yがaffine：actual F（fの三次斉次化）が被約なら
  Pε=(1+εt)∂t∧Y+Jac(f)を使う。
  F既約ではx0=1+εtの2H+cubic模型のB=εY1はtraceε≠0なので全nonzero-SS定理、
  F=LQ又は三distinct Lでは全3H+Q／5H支持定理へ入る。
  Fが非被約なら別稿のactual nonreduced cubic全source theoremを使う。
  よって**全affine Y source**はknown有限閉和Kへ入る。f3被約とは仮定しない。
* triangular Y=(A2(q)+(1−d)r)∂r+U_aff(q)、f=r g_aff(q)+H≤3(q)：
  全Ug=dg、UH+A2g=Hから別稿の明示的raw-long曲線を使う。
  曲線のPf=(r−ε²t²/2)g+Hは次数≤3であり、marked hは全εで重複したまま。
  元MULTIPLE-Hのwhole17 long／short全既知被覆の仮定を満たす。
  simple-Hのrank一modularだけからwhole17へ送る議論ではない。

## 2. 全quadratic syzygyの三係数quotient

f3=x²yの全quadratic syzygyは2yY2x+xY2y=0から

```
Y2x=xL, Y2y=−2yL,
Y2r=3b*y*r+(c/2)r²+Q2(x,y),
L=ax+by+cr.
```

と書ける。実際div零は∂rY2r=3by+crである。
linear HamiltonianによりaとQ2のx²/xy二係数を消せる。
これは全fによるactual t-shearであり、f自体を変えない。
したがって全残二次部は正確に三係数b,c,kの

```
Y2=(by+cr)(x∂x−2y∂y)
    +(3byr+c r²/2+k y²)∂r.                         (1)
```

となる。全最高syzygyを単一方向と仮定しない。

## 3. f2の全補集合

f2=A x²+Bxy+Cxr+D y²+Eyr+G r²とする。
Y1f3+Y2f2=f3のmod xは、(1)の(y,r)部がD y²+Eyr+G r²を固定する条件であり、
全四係数は

```
cG=0, 4bG−cE=0,
−4cD+bE+2kG=0, −4bD+kE=0.                       (2)
```

* G≠0では(2)がc=b=k=0を強制し、Y affineの既知枝。
* G=0,E≠0でも順にc=b=k=0となる。
* G=E=0,D≠0ではb=c=0。Y2=k y²∂rのみとなり、
  次数三のr係数からY1x,Y1yのr係数は零。
  f2にはr²/yrがなくf=r g_aff(x,y)+H(x,y)、Yqはaffineとなるからtriangular既知枝へ入る。

残るのはD=E=G=0、すなわちf2=x(Ax+By+Cr)。

C≠0ではactual r'=C r+A x+B y、t'=C t、tensor1/Cでf2=xr'へ移せる。
逆はr=(r'−A x−B y)/C、t=t'/Cであり、potentialとκ=1は保たれる。
再びlinear Hamiltonian gaugeで(1)へ戻してよい（f2は変わらない）。
このf2=xrに対する次数三比較はc=0、
Y1x=a x−k y/2−2b r、Y1y=(1−2a)yを与える。
次数二のr²係数は−2b、次のyr係数は−k/2なのでb=k=0。
従ってこの全C-unitもY affine枝である。

C=0ではactual x/y translationでf2を零にする。
x'=x+B/2、y'=y+Aの逆を元fへ入れると二次係数は全て消える。
移ったf1/f0とY1/Y0の全係数を任意として次節で処理する。

## 4. f2=0の全低次方程式

f=x²y+p x+q y+s r+d0とする。
次数三／divから

```
Y1x=a x, Y1y=(1−2a)y, Y1r=l x+m y+a r,
Y0x=−p b/2, Y0y=0, Y0r=w0.
```

全残恒等式は正確に

```
p c=0, c s=0, −2b q+k s=0, −2c q+3b s=0,
(a−1)p+l s=0, −2a q+m s=0, (a−1)s=0,
d0=−p²b/2+w0 s.                                  (3)
```

これはYf=fの次数二、一、零の全係数であり、逆に(1),(3)から全式が成立する。

* s≠0ならc=b=k=0、Y affine。
* s=0,pq≠0はa=1かつa=0を強制し空。
* s=p=0,q≠0ではb=c=0でY2=k y²∂r、fもr非依存なのでtriangular既知枝。
* s=p=q=0ではd0=0、actual F非被約の既知枝。
* s=q=0,p≠0、d0=0ではb=c=0でtriangular既知枝。

したがって最後の全残余はs=q=0、p d0≠0のみである。
(3)はa=1、c=0、b=−2d0/p²を与え、全sourceは

```
f=x²y+p x+d0,
Yx=bxy+x+d0/p, Yy=−2by²−y,
Yr=(1+3by)r+k y²+l x+m y+w0.                       (4)
```

k,l,m,w0の全四係数を保持する。新しいcompatibilityは残らない。

## 5. p,d0と四磁気係数のactual正常化

p d0≠0でa0=p/d0、b0=d0/p²とし

```
x'=a0*x, y'=b0*y, r'=a0*r, t'=t, tensor1.
```

det=1/d0、逆は非零対角の逆である。
Jac potentialはf'=x'²y'+x'+1、Yqは

```
Yx=1+x−2xy, Yy=4y²−y,
Yr=(1−6y)r+k y²+l x+m y+w0                         (5)
```

へ移る。新四係数はそれぞれ
k'=p⁵k/d0³、l'=l、m'=p³m/d0²、w0'=p w0/d0。
全四空間の可逆写像である。

任意s0≠0、A0,B0,C0について次のactual変換を取る：

```
r'=s0*r+A0*x+B0*y+C0,
t'=s0*t−2A0*y, x'=x, y'=y, tensor scalar1/s0.
```

det=s0²、逆r=(r'−A0x−B0y−C0)/s0、t=(t'+2A0y)/s0。
全六括弧への代入により(5)を保ち、係数は厳密に

```
k'=s0*k+10B0, l'=s0*l,
m'=s0*m−2B0+6C0, w0'=s0*w0−A0−C0.
```

B0=−s0k/10、C0=−s0(m+k/5)/6、A0=s0w0−C0でk',m',w0'を全て零に出来る。
l≠0では更にs0=1/lとしてl'=1。l=0は同じ全parameter多項式族のl≠0稠密開の閉境界なので、
その像も同じactual orbitの閉包へ入る。

## 6. 唯一の未帰属canonical orbit

最後に残る元は座標(t,x,y,r)で

```
Ptx=1+x−2xy, Pty=4y²−y, Ptr=(1−6y)r+x,
Pxy=0, Pxr=−x², Pyr=1+2xy.                          (6)
```

signed Pf=x²y+x+1、divP=−∂t。
actual D=h²(x²y+xh²+h³)、cubicは絶対既約（yについてprimitive一次）で、
plane cubicはcusp、h∞はcuspを通るsecantである。
基点原点の定数部はJ=∂t∧∂x+∂y∧∂rで非退化。

全元double-line最高項sourceはK∪closure(PGL5·[(6)])へ入る。
このcanonicalの既知親／新成分性は未証明。
一つの有限actual orbitへの還元であり、元sourceを一つの任意witnessで代用したのではない。
今後の有限判定は(6)一つを事前固定して行えば十分で、別点／profile／高次jetを追加する理由にはならない。
triple-line最高項x³、unimodular A=c=0全h⁵は別の残余である。
