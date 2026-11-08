# triangle double h=L^2M：二fixed源は空、残parameterは有限壁へ

2026-11-07。主義務G/E。rank二Poisson plane、f=xyz、normal位数二、
constant double G2=unit*x^2y*u^2、actualD既約被約を扱う。
全leading削減は `PLANE_M2_DOUBLE_ROOT_TRIANGLE_LEADING_REDUCTION.md` 訂正版を使う。
本稿は新しいmatrix/profile/GB/jetを使わないpureproofで、独立読解前として保存する。

## 1. 全三leading源とcanonical strict u-derivation

Delta=X-Y=diag(-1,2,-1)、Z=x partial_z、normal bracket kappa uv、k20=0。
全P1のbad u-rowはc v^2 partial_i、c nonzero。
全三leading sourceは

| bad direction | 全Y | kappa |
|---|---|---|
| x | diag(-2r-1,r+1,r) | 3r+1、r not equal -2 |
| y | diag(0,2,-2) | -2 |
| z | diag(-2,3,-1) | 1 |

であり、全raw36/忠実30 P1及び全九P2を保持する。
Rをrow divergence、Pi=P+(E wedge R)/5をcanonical strict liftとする。
全支持からR_u=kappa u、R_v=-kappa v、R_base=X+Y+span(u,v)。
従ってXi=(Pi sharp du)/uはuを固定し、vに-3kappa v/5を掛け、
base coordinatesに対してはlinear field plus bad c v^2/uだけを持つ。
Xi(F)=0はcanonical strict Casimir恒等式である。

5Xi/3のbase diagonal weightsは

\[
 \mu_j=s_j+\tfrac43d_j+\tfrac13\kappa,
 \qquad \nu_v=-\kappa,\quad\nu_u=0,
\]

ここでsjはY、djはDeltaのdiagonal entries。
unknown linear off-blockはu/vからbase座標への項だけである。
base weightsが0及びnu_vと異なる場合はactual base-from-normal shearsで全て消せる。
bad termは対応base weightが2nu_vならsemisimple部とcommuteし、
対応しない場合はdense u nonzero上のrational quadratic shiftで除ける。
その後Fをpolynomial ringへ戻すときはu-adic条件を必ず使う。

## 2. fixed bad-y源は全空

weights(x,y,z,v;u)=(-2,4,-4,2;0)。
全base weightは0及び2と異なり、normalへのlinear off-blockはactual線形shearで消せる。
bad-y termだけがresonantで、derivationは

\[
 S+N,\quad S=-2x\partial_x+4y\partial_y-4z\partial_z+2v\partial_v,
 \quad N=cv^2/u\,\partial_y,
\]

c nonzeroとなる。SとNはcommute、NはLaurent polynomial ring上locally nilpotent。
weight分解から(S+N)F=0はSF=NF=0を与える。
actual linear shearsはSplaneを保つのでFの最初のnormal項はunit*x^2 y u^2。
Nはnormal orderを一増やし、そのF2への作用はnonzero* x^2 u v^2。
F3以降のN像はnormal order四以上なのでこれを消せない。
従って全源が空。

## 3. fixed bad-z源はessential四の一actualDへ

weights(x,y,z,v;u)=(-3,6,-2,-1;0)。
全base weightは0及び-1と異なるのでlinear off-blockは同様にactualshearで消える。
resonant nilpotent部はN=cv^2/u partial_zだけ。
NF=0からpartial_z F=0で、Fはx,y,u,vだけに依存する。

S-weight零かつdegree五の全polynomial monomialsは
u^5、u^2x^2y、xyv^3だけである。
normal leadingによりu^2x^2y係数はnonzero。
xyv^3係数零ならu因子、u^5係数零ならxy因子が生じるため、irreducible部分では全三係数がnonzero。
actual diagonal GLとscaleで

\[
 F_*=xy(u^2x+v^3)+u^5\tag{1}
\]

に戻る。係数をA,B,Cとすればold x=X、old u=U、old y=(C/A)Y、
old v=tV、t^3=A/BでF=C F*、全inverseが存在する。

## 4. F*の全低degree syzygiesは一field span

base四変数を(x,y,u,v)、H=x(xu^2+v^3)とする。
degree d<=2のgradient syzygy(a,b,c,d0)を取る。
y零でHb+5u^4 c=0、gcd(H,u)=1かつdeg H=4より、c,bは共にyで割れる。
c=yL、b=yMと書くと、x零でv^3 a+5u^4 L=0。
deg L<=1なのでLはxで割れ、aもxで割れる。
L=ell*x、a=xQ（ellはdegree二の場合だけconstant）を代入し、次にx=v=0で比較するとell=0。
残式は

\[
 xu^2(2Q+M)+v^3(Q+M)+3v^2d_0=0.
\]

v零の比較から2Q+M=alpha*v。
vで割って再びv零へ制限すればalpha=0、よってM=-2Q、d0=vQ/3。
従ってAR0=0、AR1=CT、AR2=S1*Tで

\[
 T=(3x,-6y,0,v).
\]

P4のvertex変数zを戻して任意quadratic bivector Pi sharp dF*=0を考える。
各base-four rowをzの次数ごとに展開すると、係数は上のAR0/1/2に入る。
従って全imageのbase projectionはpointwise Tの一dimensional spanに入り、vertex方向を加えても二以下。
Piのrankは高々二であり、actualDがF*のrank四Poissonにはなれない。
これでfixed bad-z源も全空。

## 5. oneparameter源のgeneric whole部分

bad-x源のweights(x,y,z,v;u)は

\[
 (-r-2,\;2r+4,\;2r-1,\;-3r-1;\;0).
\]

r not in {-2,-1,0,1/2}なら全linear off-blockをactualshearで除去出来る。
bad poleはr nonzeroでrational x'=x+lambda v^2/uにより除去出来、lambda nonzero。
Fはx'のdegree高々五、uの最小Laurent exponent-5のhomogeneous degree五formとなる。

有限wall集合Wを次のように明示的に定義する。
上記四例外値に加え、nonnegative integers a,m,n,bで
a+m+n+b<=10、a+m+n<=5を満たす全tupleについて

\[
 w_r=(-a+2m+2n-3b)r+(-2a+4m-n-b)
\]

がnonzero polynomialである場合の全zeroを加える。
これは有限rational集合で、全possible Laurent monomial weightsを覆う。
根が存在しないconstant polynomialは壁を加えない。
r outside Wではweight-zero monomialの条件は恒等的weight零と同じになり、

\[
 b=n,\qquad a=2m-n
\]

である。対応Laurent monomialはu^(5-3m-n) x'^a y^m z^n v^n。
固定y/z exponentsごとに一つしかなく、rational x'を元xへ戻した最小u exponentは5-5m。
lambda nonzeroなのでm>=2はpolynomial性に反する。従ってm<=1だけ。
全Fは

\[
 A u^5+y[B(ux+\lambda v^2)^2+C(ux+\lambda v^2)zv+D z^2v^2].
\]

normal leading unit*x^2y*u^2によりB nonzero、C=D=0。
A=0なら可約なので、actuallinear scalesで

\[
 F_{\rm gen}=y(ux+v^2)^2+u^5.\tag{2}
\]

whole r outside Wの全actualDがこの一型になる。
generic source parameter又はfinite samplesを全rへ外挿しない。

## 6. Fgenもrank四を出さない

Q=ux+v^2、base四変数(x,y,u,v)とする。
degree<=2のgradient syzygy(a,b,c,d0)はQ|cを満たす。
degree二ではc=kQと書け、y零の比較Qb+5k u^4=0からk=0。
lower degreeでもc=0。従ってb=yLで、残式は

\[
 2u a+4v d_0+Q L=0.
\]

coprime(u,v)より、あるdegree<=1のMに対し
2a+xL=vM、4d0+vL=-uM。
従ってAR0=0、AR1とAR2は全て

\[
 T=(2x,-4y,0,v),\qquad Z=(2v,0,0,-u)
\]

のpolynomial multiplesに尽くされる。
vertexzを戻したquadratic Casimir bivectorのimageはbase二span plus vertex一で三以下、
skew rankは高々二。Fgenのrank四actualD源は空。

## 7. 完了範囲と明示残余

二fixed h=L^2M源は全空。
oneparameter源の全r outside有限明示Wも空。
更にr=0ではweightresonant bad-x nil部cv^2/u partial_xがF2=x^2y u^2を
nonzero normalorder三へ送り、他のresonant linear部はz方向だけなのでF2には作用しない。
第2節と同じleading-order議論でこの壁も空。

残るr in W minus {0,-2}の全wholeparameter supportは未完了。
r=-2はnormal leading amplitude零でm三以上の別義務。
壁にあるgeneric sourceを選ぶだけで全wallを覆ったとは呼ばない。
新profilefarmへ進む前に、各wallのactual-sem-plus-nil invariants又は別全source幾何を使う必要がある。
本稿は数学的還元で、まだnormal-double全triangle G完了ではない。
