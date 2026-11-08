# constant modular：最高triple-lineの全sourceを一actual orbitへ有限還元

2026-11-08、PRIMARY `E_multiple_hyperplane`、担当 `/root/plane_collision_review`。
**純全係数・正常化証明。新算術/operator/点/minor/profile/jet/GB無し。**
A=0,c≠0、最高f3がtriple lineの全源を、既知閉和Kと一つの明示actual orbitへ還元する。
最後のcanonical orbitの既知親／成分性は未解決。有限化をE所属と称さない。
producerの独立受理は別途必要である。

## 1. 全symbolと使う既知枝

κ=1、base順序(x,y,r)でP=∂t∧Y+Jac(f)、divY=1、Yf=f、degY≤2を用いる。
f3のunused base方向の尺度を使ってf3=x³にactual正規化し、全低次を残す。
Y2x=0、divY2=0からplanar polynomial Poincareにより

```
Y2=H3_r∂y−H3_y∂r,
H3=x Q2(y,r)+Q3(y,r)                               (1)
```

と出来る。元のH3のx²y/x²r項は全fによるlinear Hamiltonian t-shearで消し、x³は場を変えない。
Q2は全三係数、Q3は全四係数であり、全七係数quotientを保持した。
既知枝はactual fがcubeとなる全源と、全triangular sourceである。
後者は別稿のraw-long actual曲線で、全Pf次数≤3／marked h重複を保持したwhole17 theoremへ入る。

## 2. Q3がcubeでない非零式ならactual fはcube

homogeneous binary cubic Aとquadratic BについてdA∧dB=0なら
Euler縮約から3A dB−2B dA=0、したがってB³=cA²。
B≠0ならUFDの各valuationからAはlinear cubeでなければならない。
linear first integralも同様にA=cL³を強制する。

次数三のmod xから、Q3のHamiltonian場はf2|x=0を固定する。
Q3がcubeでないならf2=xL1。
次のmod x²からQ3のHamiltonian場はL1の(y,r)部を固定するのでL1=a x。
actual x translationでf2=a x²を消す。
今度は次数二のmod xからf1の(y,r)部も零。
残るY1x=x/3、Y0x=0とYf=fの次数一・零からf1=f0=0。
よってactual f=x³で、全sourceは既知cube枝である。
この議論はQ3のdouble+simple型／三単根型と全Q2に同時適用する。

## 3. Q3=0：Q2の全三rank

actual GL2(y,r)、非零尺度、linear Hamiltonian gaugeでQ2の全rankを扱う。
q→Aq+xvのshearが(1)に加えるx²q項も同じgaugeで消せる。

### rank二：H3=xyr

Y2=x(y∂y−r∂r)。
f2=A x²+x(By+Cr)+D y²+Eyr+G r²と置くと、次数三からD=G=0、
Y1x=x/3−By/3+Cr/3。
E=0なら次数二のy²/r²係数が−B²/3、C²/3なのでB=C=0。
x translationで残るAx²を消すと、次数二からf1のy/r係数とY0xが零、
次数一・零からactual f=x³となる。

E≠0ならq→q+xvでf2のxy/xrを消し、x translationでx²を消してf2=Eyrにする。
最高Y2は上記gaugeで元へ戻せ、次数三からY1x=x/3。
divY1=1より(y,r)blockのtraceは2/3である。
次数二のyr係数はY1f2だけから2E/3、一方右辺はE。
Y2f1とY0f3はxを含むので寄与せず、E≠0に反する。全E-unitは空。

### rank一：H3=xy²

Y2=−2xy∂r。次数三からE=G=0、Y1x=x/3+(2C/3)y。
C≠0,D=0では次数二のyr係数が2C²/3で、他項から打ち消せず空。
C D≠0ではf2の対称matrixは

```
S=[[A,B/2,C/2],[B/2,D,0],[C/2,0,0]], detS=−D*C²/4.
```

S⁻¹のxx/xy成分は零。次数二の全式Y1f2+Y2f1+Y0f3=f2にS⁻¹ traceを取ると、
Y2f1はxy倍、Y0f3はx²倍なのでtrace寄与零。
従って2trY1=3、divY1=1と矛盾する。全C-unitが空である。

C=0ならf2=A x²+Bxy+D y²。
B又はDが非零なら、次数二のr係数からY1yのr係数が零。
Y2はvertical、fはrについてaffineなので全triangular sourceへ入る。
B=D=0ならAx²を消し、次数二からf1のr係数が零。
f1のy係数非零なら次数一から同じtriangular構造、y係数零ならactual f=x³。

### rank零

Y2=0なので全affine-Y既知枝を使う（actual F被約／非被約を区別する既補題）。

## 4. Q3がcube：全Q2の二係数正常形

Q3≠0がcubeならactual frame/尺度でQ3=y³と出来る。
y→y+αxによりQ2のy²係数を消し、現れるx²q項をHamiltonian gaugeで消す。
全源は

```
H3=y³+x(byr+c r²),
Y2y=bxy+2c x r, Y2r=−3y²−b x r.                   (2)
```

b,cの全二係数を保持する。
f2=A x²+Bxy+Cxr+D y²+Eyr+G r²に対し、mod xでE=G=0、
続く次数三から

```
C=2Db/3, D c=0,
Y1x=x/3−Bb*y/3−(2Bc−Cb)r/3.                       (3)
```

c≠0ではD=C=0。
Ax²をtranslationで消した後、次数二のyr係数−2B²c/3からB=0。
残る全次数二・一・零からactual f=x³となる。

c=0,b=0ならY2=−3y²∂r、C=0。
B又はDが非零なら次数二r係数からY1yのr係数が零でtriangular源。
B=D=0ならAx²を消し、f1のr係数は次数二から零。
f1のy係数非零では次数一からtriangular、残りはactual cubeである。

c=0,b≠0の場合、D≠0ならC=2Db/3≠0であり、
次数二のr²係数C²b/3が非零となって矛盾。
従ってD=C=0。Ax²を消し、B=0の場合も全低次からactual cubeへ入る。
以下でbB≠0の全残余を処理する。

## 5. bB非零の全残source

f=x³+Bxy+p x+q y+s r+dと書く。
全次数二／一／零の比較は、s=−B²b/9、q=B/(2b)を強制し、
他の全自由係数を次節で明示するsourceへ送る。
この導出は次数二から

```
Y1x=x/3−Bb*y/3,
Y1y=l x+m y+n r, Y1r=L x+M y+N r,
n=bs/B, m=2/3−bq/B, N=bq/B,
Y0x=−B*l/3,
s=−B²b/9
```

を得て、次数一のr係数q n+s N=sからq=B/(2b)を得るものである。
s≠0なので除算は元のbB-unitに従う。
残る次数一のx/yと次数零はY0y,M,Y0rを一意に定める。

a=−b、c0=b²B/3、s0=c0³/a²としactual尺度

```
x'=a x, y'=c0 y, r'=c0² r, t'=s0 t, tensor1/s0
```

を取る。全係数は非零で両逆は対角逆、
b'=−1、B'=3、q'=−3/2、s'=1となり、f3=x³と(2)のy³係数も保つ。
p'=a²p、d'=a³d。残る全linear係数もactual push後に上の全解へ戻る。

以後記号を戻し、全残源はp,d,l,Lの四変数による必要十分表示

```
f=x³+3xy+p x−3y/2+r+d,
Yx=x/3+y−l,
Yy=−xy+l x+y/6−r/3+v,
Yr=−3y²+xr+L x+(3l−p−5/4)y+r/2+w,
v=2p/9+l/2−L/3,
w=d+lp+p/3+3l/4−L/2.                              (4)
```

となる。代入によりdivY=1,Yf=fは全四係数で恒等的に成立し、
逆も上の全係数比較で得た。未知compatibilityやthin境界を捨てていない。

## 6. 全四係数を消すactual逆付き正常化

まず

```
r'=r+p x+d, t'=t+(p/3)y, x'=x,y'=y, tensor1
```

とする。det1、逆r=r'−p x−d、t=t'−p y/3。
全fはf0=x³+3xy−3y/2+rへ移り、同じ(4)で
p'=d'=0、l'=l+p/3、L'=L−d−p/6となる。

残るl,Lを再びl,Lと書く。
任意Aに対し次のactual affine f0-automorphismを取る：

```
X=x+A, Y=y−A x,
R=r−(3A/2)x−3A y−A³.
```

f0(X,Y,R)=f0(x,y,r)、det1。
逆は

```
x=X−A,
y=Y+A X−A²,
r=R+(3A²+3A/2)X+3A Y−2A³−3A²/2.
```

次にtを

```
T=t+βY+γR+αX,
β=−A²−A/2, γ=−A/3
```

へshearする。逆はt=T−βY−γR−αX。
base変更で最高Hamiltonian potentialに現れる
−A x²r−(3A²+3A/2)x²yはγ/βの全f-Hamiltonianで正確に消える。
全係数比較により(4)のp=d=0を保ち、

```
l'=l−2A/3,
L'=L−2A³−2A²−A−3α.
```

従ってA=3l/2、α=(L−2A³−2A²−A)/3で全l,Lを零に出来る。
どのl/Lも反転せず、全四係数域が一つのactual orbitへ移る。
全変更はaffine/PGL、tensor尺度は先の非零対角尺度だけで、非線形Darbouxは用いない。

## 7. 唯一残るcanonicalと正確なscope

最後の元は

```
f0=x³+3xy−3y/2+r,
Y=(x/3+y)∂x+(-xy+y/6−r/3)∂y
  +(-3y²+xr−5y/4+r/2)∂r,
P=∂t∧Y+Jac(f0).
```

全六括弧は

```
Ptx=x/3+y, Pty=−xy+y/6−r/3,
Ptr=−3y²+xr−5y/4+r/2,
Pxy=1, Pxr=−3x+3/2, Pyr=3x²+3y.
```

divY=1、Yf0=f0、Pf=f0を直接手展開出来る。
actual D=h²[x³+3xyh−(3/2)yh²+r h²]。
cubicはrについてprimitive一次式で絶対既約、rank四である。
必要なら同じtensorのactual affine基点r=1でPf=1を取れるが、新しいsource点を選ぶ理由にはしない。

全triple-line最高項sourceは、既知有限閉和Kとこの**一つのactual orbit**で覆われる。
親への所属／成分性は未証明。double-line最高項の別canonicalとはactual cubicのessential変数数が異なり、
同点であるとの推測はしない。今後はこの元一つを固定した最弱のknown-parent判定で十分である。
unimodular A=c=0の全h⁵源は本稿の外。
