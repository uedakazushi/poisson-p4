# constant modular：irreducible nodal/cuspidal最高三次式の全MPS被覆

2026-11-08、PRIMARY `E_multiple_hyperplane`、担当 `/root/plane_collision_review`。
**純係数証明。新operator・点・rank・minor・kernel計算・jet・GBは実行しない。**
A=0,c≠0の全源のうち、fの最高三次式f3が既約かつ特異な全支持を処理する。
f3がreducible/nonreducedの補集合は本稿では未処理として保持する。
独立pure査読/root採録はproducerの本稿とは別である。

## 1. 全低次項を保つ入口

`MULTIPLE_HYPERPLANE_LOW_RANK_MODULAR.md` §3 の必要十分表示に従い、
非零tensor scalarとactual affine frameでκ=1として

```
P=∂t∧Y+Jac(f), divY=1, Y(f)=f, degY≤2, degf≤3.
```

Y=Y2+Y1+Y0、f=f3+f2+f1+f0と全ての低次項を残す。
次数四とdivの次数一から

```
Y2(f3)=0, divY2=0.
```

を得る。以下でこの全二条件を手で解き、Y2が線形Hamiltonianに尽くされることを示す。
一般的なsmooth f3のKoszul完全性をsingular f3へ外挿しない。

## 2. 既約特異三次式の全二型

nodal cubicはnodeを(0:0:1)、二接線をx=0,y=0にして
F=xyz+A x³+B x²y+Cxy²+D y³と書ける。
zのlinear shearでB,Cを消す。既約性からAD≠0で、可逆尺度により
F_N=x³+y³+xyzへ移る。
cuspidal cubicは二重接線をx=0としてF=x²z+binary cubic(x,y)。
zのlinear shearでx³,x²yを消し、yのshearと再びzのshearによりxy²を消せる。
y³係数が零なら可約なので、可逆尺度によりF_C=x²z−y³へ移る。

したがって全既約特異f3はactual GL3変更後にcF_N又はcF_C、c≠0となる。
残るcを一へ固定したとは仮定しない。syzygy空間とtrace条件はcで変わらず、
Hamiltonianを戻す際だけ対応する線形関数を1/c倍すればよい。

## 3. cuspの全quadratic syzygyとdivergence

F=x²z−y³とし、quadratic vectorを(A,B,C)とする。
全方程式は2xzA−3y²B+x²C=0である。
mod xからB=xL、もう一度mod xとgcd(z,y)=1から
L=mx+l z、A=(3l/2)y²+x(nx+py+qz)を得る。
残るCは一意で、全解は

```
A=3l*y²/2+x*(n*x+p*y+q*z),
B=m*x²+l*x*z,
C=3m*y²−2z*(n*x+p*y+q*z).
```

すなわちm,l,n,p,qの全五係数であり、追加のchartやunitを使わない。
div(A,B,C)=−py−3qzなので、div零は正確にp=q=0。
このとき

```
Y2=∇F×∇(m*x−n*y−l*z/2).
```

全divergence-free quadratic syzygyは三つの線形Hamiltonianだけとなる。
同じ係数比較をlinear vectorへ適用すると、Fを保存する全linear場は
n(x∂x−2z∂z)だけで、そのtraceは−nである。
したがってtrace零でFを保存するlinear場は零。

## 4. nodeの全quadratic syzygyとdivergence

F=x³+y³+xyzとする。
全18係数のquadratic vector(A,B,C)へ∇F·(A,B,C)=0を課し、
x⁴,y⁴,xz³,yz³,x²z²,y²z²から順に比較すると、全解は次の四係数b,c,h,kである：

```
A=b*x*y+c*x*z−3k*y²,
B=−3c*x²+h*x*y+k*y*z,
C=−3b*x²+9(c+k)*x*y−h*x*z−3h*y²−b*y*z−(c+k)*z².
```

逆にこの式を∇Fへ掛けると全quartic係数が相殺するため、必要十分な全表示である。
特にdiv(A,B,C)=−(c+k)z。
div零はk=−cと同値で、その全解は

```
Y2=∇F×∇(h*x−b*y+c*z).
```

となる。三つのHamiltonian以外のsyzygyを捨てたのではなく、そのdivergenceを課して除いた。
linear vectorの場合はx³,y³,xz²,yz²等を比較すると全係数が零となり、
Fを保存する非零linear場は存在しない。

## 5. 全低次項のactual消去

上の各全Y2を、必要なc尺度も含めY2=∇f3×∇ℓと書く。
ℓはlinearである。actual affine変更T=t−ℓ(x,y,z)はdet1、逆t=T+ℓ。
新YはY−∇f×∇ℓなので、最高二次部が全て消える。
Hamiltonian場はdivergence零でfを固定するから、divY=1、Yf=fは全低次項込みで保たれる。
非線形Darboux変換をPGLとして用いていない。

以後Y=Y1+Y0。次数三からY1f3=f3、divY1=1。
Y1−Euler3/3はf3を保存しtrace零なので、§§3–4により零である。
従ってY=Euler3/3+Y0。
actual translation x'=x+3Y0によりY=Euler3'/3となり、同じ全多項式fはYf=fから
純三次式になる。f2/f1/f0は初めから捨てたのではなく、この恒等式により全て零となる。
よって元の全源はactual PGLで

```
P=∂t∧Euler3/3+Jac(g3)
```

へ移り、g3は全既約node/cusp cubicを含む。

## 6. 全MPS閉和への実族

任意のternary homogeneous cubic gとεに対して

```
P(g,ε)=(1+εt)∂t∧Euler3/3+Jac(g)
```

は真の二次Poisson族である。
Euler3(g)=3g、divEuler3=3からEuler3がJac(g)を保存し、t方向の追加も全Jacobiを満たす。
Pf=(1+εt)g、actual五次D=h(h+εt)g。
g=xyz、ε≠0では五つの独立hyperplaneによるSNC点となる。
全(g,ε)のparameter空間は既約で、その像閉包を含む全Poisson成分はこのSNC点を含む。
既受理MPSのSNC-meeting component定理によりその成分はMPS40の一つである。
従って全族、特にε=0の全gはMPS有限閉和へ入る。

§5の可逆frame/shear/translationを戻すことで、全既約nodal/cuspidal f3源を被覆する。
新成分増分零。全sourceの単位条件はc≠0と元κ≠0だけで、その条件はf3既約三次／
constant modular非零という入口そのものである。全未知低次fiber、primitive β、
smoothなf3の一般点からのsource平坦性は仮定しない。

## 7. 保持する残余と反例

このdivergence-free syzygy結論は全reduced cubicへ無条件に拡張出来ない。
例えばf3=xy(x+y)は三concurrent linesで、任意binary quadratic Q(x,y)∂zが
syzygyかつdivergence零であり、全三係数のうち∇f3による二次元Hamiltonian像以外も残る。
本稿はこのconcurrent枝、line+conic、他のreducible/nonreduced f3を処理したと主張しない。
またA=c=0のunimodular全h⁵ sourceも本稿の外である。
