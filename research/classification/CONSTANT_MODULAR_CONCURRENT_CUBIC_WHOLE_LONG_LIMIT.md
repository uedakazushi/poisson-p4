# constant modular：concurrent三線最高項の全sourceとactual long極限

2026-11-08、PRIMARY `E_multiple_hyperplane`、担当 `/root/plane_collision_review`。
**純証明。新operator・点・minor・kernel計算・profile・jet・GB無し。**
A=0,c≠0、f3が相異なる三concurrent linesである全低次sourceを扱う。
whole longの特殊な非線形straighteningをPGLと呼ぶことなく、元六括弧への
actual affine曲線から帰属を証明する。独立査読/root採録は別途必要である。

## 1. 全最高syzygyと余剰一方向

元必要十分入口をκ=1、座標(t,r,a,b)で

```
P=∂t∧Y+Jac(f), divY=1, Yf=f, degY≤2, degf=3
```

とする。Jacの三変数順序は(r,a,b)。
f3はactual linear frameにより非零係数cを保持して
F(a,b)=c*ab(a+b)となる。Fa,Fbは互いに素なbinary quadraticsであり、
Y2=(A,B,C)を(a,b,r)順に書けばAf_a+Bf_b=0から
(A,B)=λ(Fb,−Fa)、λ定数である。
この二成分のdivは零なのでdivY2=0はCのr非依存性と同値。
したがってCは任意binary quadraticの全三係数である。

linear Hamiltonianは、r方向のlinear関数で上記λを消し、
(a,b)方向のlinear関数でverticalのspan(Fa,Fb)二次元を消せる。
actual t-shear後に残る全二次部は

```
Y2=μ Q2(a,b)∂r,    Q2 not in span(Fa,Fb)
```

だけである。Fa,Fb,Q2はSym²(a,b)の基底。例えばQ2=abを選べる。
μ=0はYがaffineとなるので、別稿のaffine-Y実pencil補題が全sourceへ適用される。
以後μ≠0を保持し、その補集合は既処理のμ=0である。

## 2. μ非零の全低次係数がtriangular sourceへ入る

一般のY1のa,b成分におけるr係数をp,qとし、f2のr²係数をsとする。
Yf=fの次数三のr係数は厳密に

```
p Fa+q Fb+2s μQ2=0.
```

三つのbinary quadricが基底でμ≠0なのでp=q=s=0。
この全零条件により、他の係数を特別化せず

```
Y=[A2(a,b)+(1−d)r]∂r+U(a,b),   d=div U,
f=r g(a,b)+H(a,b)
```

となる。Uは全affine場、A2は任意二次以下、gはaffine、Hは三次以下でH3=F。
上の全低次項からfのr²が消え、Yのa,b成分もr非依存となった。
Yf=fの全係数比較は正確に次の二恒等式である：

```
U(g)=d*g,
U(H)+A2*g=H.
```

divY=1は(1−d)+d=1で恒等的に成立する。
この表示は元の全Pを含む必要sourceであり、逆にこの二式から元Jacobi/div/Pf条件も成立する。
したがって未知compatibilityや低次元sourceを途中で捨てていない。

## 3. 任意triangular sourceからraw longへのactual曲線

ここでは§2の二恒等式だけを仮定する。旧raw long source(19)の座標を(X,Y,a,b)とし、
ε≠0で次を代入する：

```
Uold=−U, dold=−d,
Fold=A2, gold=g/ε, Hold=H/ε.
```

全条件Uold(gold)=dold*gold、(Uold+1)Hold=Fold*goldは
上の二恒等式から厳密に従う。次数上界も全て保持される。
旧sourceのr_old=Y−X²/2は元式を記述するためだけの略記である。
次の変更はそれとは異なり、actual affine変更である：

```
t=−X/ε, r=Y, a=a, b=b; tensor scalar ε, Euler0.
```

detは−1/ε、逆X=−εt,Y=rである。
全六括弧は元raw式から直接

```
Ptr=A2+(1−d)r−(1−d)ε²t²/2,
Pta=Ua, Ptb=Ub,
Pra=ε²t*Ua+(r−ε²t²/2)*g_b+H_b,
Prb=ε²t*Ub−(r−ε²t²/2)*g_a−H_a,
Pab=g.
```

となる。全係数がε多項式で次数≤2、ε=0は§2の元tensorそのものである。
全Jacobiはε≠0のactual long所属から成り、同じ多項式恒等式として全εへ延長する。
signed Pfは正確に

```
fε=(r−ε²t²/2)*g+H.
```

これは六括弧の直接Pfでも、oldPf=−(r_old*gold+Hold)、
座標det−1/εとtensor scalar²=ε²のcovarianceからも得られる。
H3=c*ab(a+b)≠0なのでPfは全εで非零。全rank四点を保ったactual曲線である。

## 4. 曲線の親は既成分二つだけで十分

gが非定数なら、ε≠0でgoldも非定数。
旧longの非定数g全源は `RANK_ONE_LONG_NONCONSTANT_25_COMPONENT.md` の
既受理25成分に全点帰属する。全U/H/F/frameの補集合を含む定理を使う。
従ってε=0の元tensorも同じ閉成分に入る。

gが非零定数ならU(g)=dgからd=0。
一方、H3=Fがbinary squarefreeなので、U1F=FはU1=Euler2/3を強制する。
実際、Fa,Fbに一般2×2linear場を掛けてa³,b³,a²b,ab²を比較すれば、
off-diagonal二係数が零、二対角が1/3と得る。
よってd=2/3となり矛盾。この全枝は空である。

g=0ならUH=Hであり、同じ最高係数比較からU1=Euler2/3。
actual translation q'=q+3U0によってU=Euler2'/3へ移り、
UH=HからHは中心に関する純binary cubicとなる。
その三次最高項はsquarefreeなので、この全sourceは旧longの
SC：Uold=−Euler2'/3、gold=0、Hold全binary cubic、Fold全quadraticへ入る。
SCの全source閉包が25次元成分であることは
`RANK_ONE_LONG_THREE_AMBIENT_SEPARATE_REVIEW.md` と全源pure入口の独立査読により受理済み。
従って全ε≠0とε=0はその同じ閉成分に属する。

以上はg非定数／非零定数／零の全三補集合を尽くす。
新しいM1/M11成分性を本稿の前提へ追加する必要はない。

## 5. 全最高被約三次式の合成と残余

μ=0の全affine-Y源、μ≠0の全triangular源、全g三枝、全low termsと元shear/frameを戻し、
最高項が三相異concurrent linesである全sourceがknown有限閉和へ入る。
別稿のreduced essential3最高項の全被覆と合わせれば、
**A=0,c≠0でf3被約の全支持**はknown閉成分で覆われる。
この合成は両producer純証明の独立受理を条件とする。

f3がdouble-line+line又はtriple-lineの非被約の場合、本稿のgradient独立性を使わない。
その全追加syzygy、及びA=c=0のunimodular全h⁵ sourceは別の残余である。
新成分増分零、全A0の完成／全分類の完成は主張しない。
