# m2 moving double root：全weight零sourceは一つの全parameter族

2026-11-07。主義務G/E。reduced nonzero R、rank二Poisson plane、実D既約被約、
位数二のnormal quadraticがdoubleでprimitive root rowのbase degreeが一である全sourceを扱う。
`PLANE_M2_DOUBLE_ROOT_CONNECTION_GEOMETRY.md` の全接続還元を使う。
本稿はP0の完全正常形であり、未知の全P1/P2の被覆又は次元上界は主張しない。
以下は有理座標の手計算で、新しいsource/minor/GB/jet jobは無い。

## 1. 全位置と全scalar period

G2=h ell^2、ell=a u+bv、q=f/h、p=(a=b=0)とする。
既証明からp in h intersect q、hはRのline成分、qはreduced quadratic。
linear completionMのfirst rowは(a,b)、det M=qであり、全connectionは
beta=theta I-M^-1dM。thetaは全closed logarithmic scalar periodを保つ。
Mのsecond rowは補助frameで、rootのnormalizationには使わない。

## 2. qが二線なら全P0が空

q=LMとする。p=L intersect Mのとき、定数normal右GL2でell=Lu+Mvと出来る。
Mをconstant left matrixとdiag(L,M)の積に取れるためbetaはdiagonal。
actualG2のu^2係数は必ず3f k20だが、targetはhL^2のnonzero倍であり、
f=hLMで割れない。従って空。

p in L minus Mのとき、hとLはpのidealを生成するのでell=hu+Lvと出来、

\[
 M=\begin{pmatrix}h&L\\-M&0\end{pmatrix},\qquad
 M^{-1}dM=
 \begin{pmatrix}d\log M&0\\(dh-h\,d\log M)/L&d\log L\end{pmatrix}.
\]

betaはlower triangular。canonical mixed u-columnはu倍なので、actualG2のv^2係数は3f k02。
target hL^2はfで割れず再び矛盾。
L/Mを交換した補集合も同じであり、triangle及びconcurrent三線の全moving-double sourceが空になる。
marked rootをmixするleft GL2によるsource同一視は用いていない。

## 3. affine h=1 chartの正しい使い方

残るqはsmooth conic。D(P0)=hQ^2のsimple component hはP0不変なので、
h=1 chartのP0はquadratic affine Poisson tensorである。
これはP0だけの必要記述であり、全P=P0+P1+P2でh不変とは置かない。
affine base coordinatesをr,sとし、normal root frame(U,V)=M(u,v)を使う。
両rootのmixed fieldsは同じscalar Hamiltonian fieldXで、
{r,U}=U Xr、{r,V}=V Xrとなる。
det mixed blockは零なのでaffine PfとJacobian detMの変換則から

\[
 \{U,V\}=kU^2,\qquad k\ne0.\tag{1}
\]

k nonzeroはG2 nonzeroそのもの。base-zero cubicを固定して全体scaleを一度規格化する。
Rはconicとlineの二成分なのでthetaのaffine restrictionはconicのdlogのconstant倍、
従ってXは下の各一parameter全sourceである。

## 4. transverse conic+lineは全空

h=1でq=1-xz、base bracket{x,z}=3q、root U=u+zv、補助V=xu+vと出来る。
M=[[1,z],[x,1]]、det M=qであり、全scalar fieldは
X=t(x partial_x-z partial_z)。
normal bracketK={u,v}はP0 affine quadraticだから、u^2係数K20はconstantである。
Leibnizと二scalar root equationsを全て保って計算すると

\[
 \{U,V\}=qK-t x uU-t z vV+3q uv.
\]

式(1)とのu^2係数比較はq K20=k+t x。
q=1-xzなので左辺がconstant polynomialの倍であることとk nonzeroは両立しない。
従ってtransverse conic+lineも全moving-double sourceが空。

## 5. tangent conic+lineの全P0

base coordinates(h,x,y)でf=h(hx-y^2)とし、h=1 chartでq=x-y^2、
{x,y}=3q、U=yu+v、V=xu+yvと取る。
det M=y^2-x=-q、全scalar fieldはX=t(2y partial_x+partial_y)。
Leibnizから

\[
 \{U,V\}=(y^2-x)K+t(x-2y^2)u^2-2tyuv-tv^2+3qu^2.
\]

式(1)とのv^2係数及びdegree二polynomial性からk=-t、t nonzero。
uv係数は零、u^2係数はK=(t+3)u^2を与える。
全mixed bracketsを戻すと、h=1 chartの全six bracketsは正確に

\[
\begin{aligned}
 P_{xy}&=3(x-y^2),&P_{xu}&=(2t+3)yu-3v,\\
 P_{xv}&=-3xu+(2t+3)yv,&P_{yu}&=(t+3)u,\\
 P_{yv}&=-3yu+tv,&P_{uv}&=(t+3)u^2.
\end{aligned}\tag{2}
\]

Pfaffianは3t(yu+v)^2。
Jacobiはdense rational frame(x,y,U,V)で、flat scalar connectionとparallel fiber bivector
-t U^2 partial_U wedge partial_Vによるproduct constructionから全tで成立する。
同じことは表示六括弧のLeibniz展開からも直ちに確認出来る。
分母qはdense openでだけ使い、式(2)はpolynomialなのでq=0の点も恒等式に含まれる。

元のprojective P0はh-Casimir quadratic liftで得る。斉次六括弧は

\[
\begin{aligned}
 P_{xy}&=3(hx-y^2),&P_{xu}&=(2t+3)yu-3hv,\\
 P_{xv}&=-3xu+(2t+3)yv,&P_{yu}&=(t+3)hu,\\
 P_{yv}&=-3yu+thv,&P_{uv}&=(t+3)u^2,
\end{aligned}\tag{3}
\]

P_h*=0。実Dは3t h(yu+hv)^2である。
t=-3等のspecial coefficientsも除かず、全parameter domainはt nonzeroだけ。
normal rootのframeとbase cubicのscaleを戻すactual linear mapsにより全moving-root P0を覆う。

## 6. 次の必要最弱義務

式(3)は全moving-double P0の唯一の残familyである。
weights(h,x,y,u,v)=(0,2,1,N,N+1)に対しP0のcoefficient-minus-index weightは-1で、
normalweight正の全P1/P2はN sufficiently largeでrelative positiveになる。
従って未知のfull sourceはP0へのactual orbit limitを持つが、
P0の既知MPS所属だけから全higherの所属は従わない。
既 `DOUBLE_QUADRIC_SUPPORT.md` 第7節のgeneric contact-source smooth pointと比較する場合にも、
allambient unique-component仮定及びspecial-parameter補集合を別に確認する必要がある。

本稿のfullhigher義務はこの一parameter全族だけに縮んだ。
constant-root d0の全nonlog closed-form source、nonreduced R、m3及び他surfaceは別残余。
次の計算を行うなら、この全tと全boundaryに対する最弱fullhigher命題を先に契約し、
既接触sourceの無用な再計算又はpoint profile列へ進めない。
