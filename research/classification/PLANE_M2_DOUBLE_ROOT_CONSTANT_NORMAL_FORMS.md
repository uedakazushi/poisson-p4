# tangent conic及びconcurrent三線：constant double rootの全P0模型

2026-11-07。主義務G/E。実D既約被約、rank二restriction plane、reduced R、位数二、
generic double quadraticでprimitive root degree零とする。
本稿は全P0の正常形とfullhigherのbad-row必要条件だけを記録する。
未知P1/P2の被覆又はG boundは未証明であり、新計算jobは無い。
triangleは `PLANE_M2_DOUBLE_ROOT_TRIANGLE_LEADING_REDUCTION.md` 訂正済み版、
transverse conicは `PLANE_M2_DOUBLE_ROOT_CONIC_REDUCTION.md` の全空性を使う。

## 1. 全記法とactualPfを先に固定する

rootはu、G2=gamma h u^2、gamma nonzero。
raw sourceはP_f、mixed uXとuZ+vY、normal bracket k20u^2+kappa uvである。
Delta=X-Y、iE iDelta volume=-3f dlog(f/h)。
ここで扱うtangent/concurrentでは全X,Yが同じnilpotent spaceに入り、
X(f)=Y(f)=0、X wedge Y=0なのでtargetuv0からkappa=0。
全actual u^2係数は正確に

\[
 G_2/u^2=3f k_{20}-\det(E,X,Z)-Z(f).\tag{1}
\]

全Zを保持してこの式を使う。Z非零のdet項を無断で消さない。
全off-diagonal connectionはclosed omega=(f/h)b=N2/hから戻り、iE iZ volume=-3h omega。
fullhigher bad u-row v^2 CにはC nonzeroが必要である。

## 2. tangent conic：h=fは二つのactual P0軌道

f=z(xz-y^2)、N=2y partial_x+z partial_y。
既全bad-row解析でX=Y=-2N、C in C*partial_xが必要。
h=fなのでbはclosed logarithmic、全Z=cN。
式(1)ではdetとZ(f)は零なのでk20 nonzero、相対normal scaleでk20=1に出来る。

base weights(x,y,z)=(2,1,0)、normal weights(u,v)=(0,1)のactual diagonal GL5を用いる。
P_fとdiagonal mixedのweightは-1、normal u^2 termも-1、off uN partial_vのweightは-2。
projective scaleを合わせるとcだけが非零characterで変わる。
従ってc=0又は1が全P0を覆う。c=0の閉補集合を非零orbitと同一視しない。
全P1/P2はこの変換で同時に戻すが、その解集合はまだ決定していない。

## 3. tangent conic：h=z^3は一つのactual P0軌道

既bad-row解析でX=-8N、Y=-5N、C in C*partial_xが必要。
closed omegaの全五parameter（z=1 chartのdegree二potential）から全Zは

\[
 Z=(A x+B y+C_0 z)\partial_x+(D x-A y+E_0 z)\partial_y
\]

であり、z-componentは零。
式(1)を全monomialで比較すると

\[
 D=A=k_{20}=0,\qquad B=2E_0,\qquad C_0\ne0.
\]

実際det(E,X,Z)+Z(f)は

\[
 -18Dxyz+18Ay^2z+9Axz^2+(-18E_0+9B)yz^2+9C_0z^3.
\]

従ってZ=E0 N+C0 z partial_x。
normal shear v to v+c uはZをc Delta=-3cNだけ変えるためE0を消せる。
相対normal scaleによりC0を一へ戻し、全模型は

\[
 X=-8N,\quad Y=-5N,\quad Z=z\partial_x,\quad k_{20}=\kappa=0.
\]

この模型のactualG2は-9z^3u^2で非零。
従ってtangentのconstant double rootは前節の二模型と本節の一模型、計三つに尽くされる。
moving rootの別oneparameter全族は `PLANE_M2_MOVING_DOUBLE_ROOT_NORMAL_FORM.md` に残る。

## 4. concurrent三線：全bad Cは共通radical方向

f=xy(x-y)。closed diagonal fieldsは
Y=(a x+b y)partial_z、X=Y+Deltaであり、hの全三型に対し

| h | Delta=(r x+s y)partial_z |
|---|---|
| f | r=s=0 |
| x^2y | r=0,s=3 |
| x^3 | r=-3,s=6 |

全S3 marked choicesはbase線形変換で戻す。
fullweight一のvector式へC=(c,d,e)を代入すると、base componentsは次の全四式を満たす。

\[
 d(a-r+2)=0,\quad (a+2)c+(2b-s-2)d=0,
\]
\[
 (r+2-2a)c-(b+2)d=0,\quad c(b-s+2)=0.\tag{2}
\]

c nonzero,d zeroならr=-6、d nonzero,c zeroならs=-6が必要。
両方nonzeroならa=r-2、b=s-2であり、残二式の2x2 determinantは6(6-r-s)。
表の三型ではいずれも非零で、r=-6/s=-6にもならない。
従って常にc=d=0、C in C*partial_z。
これはYの全二parameterを保つ結論で、periodをgenericと仮定していない。

## 5. concurrent h=f：全四parameter leading源

全off-fieldはZ=(c x+d y)partial_z。
式(1)でdetとZ(f)は零なのでk20 nonzero。
相対normal scaleでk20=1として、全X=Y=(a x+b y)partial_z、全Z二parameterを保つ。
これ以上Aut(R)のnormalizerで消せないparameterを無理に除かない。
全P0はこの四parameter源とactual linear飽和で覆われる。

## 6. concurrent h=x^2y：全三parameter leading源

closed omegaの全三parameterから

\[
 Z=A_0 y\partial_z+B_0T+C_0 x\partial_z,\qquad
 T=\operatorname{diag}(-1,2,-1).
\]

X=(a x+(b+3)y)partial_z、Y=(a x+b y)partial_z。
式(1)は

\[
 G_2/u^2=3xy[(x-y)k_{20}+B_0a x+B_0(b+4)y].
\]

従ってB0 nonzero、a+b+4 nonzero、k20=B0(b+4)が必要。
A0とC0はこの式へ寄与せず、全係数を保持する。
normal shearでDelta=3y partial_z方向のA0を消し、相対normal scaleでB0=1に出来る。
残る全模型は

\[
 Z=T+c x\partial_z,\quad k_{20}=b+4,\quad a+b+4\ne0,
\]

a,b,cの全三parameter源である。
normal bracket k20=0のb=-4境界も許す。非零条件はk20ではなくa+b+4である。

## 7. concurrent h=x^3：全二parameter leading源

全closed omegaの五parameterから

\[
 Z=(A_0x+B_0z+3E_0y)\partial_y
 +(C_0x+D_0y-3E_0z)\partial_z.
\]

X=((a-3)x+(b+6)y)partial_z、Y=(a x+b y)partial_z。
R1=4-a、R2=-b-8と置くと、式(1)は

\[
 G_2/u^2=3xy(x-y)k_{20}
 -x(A_0x+B_0z+3E_0y)(R_1x+R_2y).
\]

target gamma x^3 nonzeroにより

\[
 A_0\ne0,\quad a\ne4,\quad B_0=0,
 \quad k_{20}=E_0(b+8),\quad
 A_0R_2+3E_0(R_1+R_2)=0.
\]

a+b+4=0なら最後の式はb=-8,a=4を強制しa not equal fourに反する。
従ってa+b+4 nonzeroであり、相対normal scaleによりA0=1として

\[
 E_0=-\frac{b+8}{3(a+b+4)}.
\]

残るC0,D0はactual base shears z to z+p x+q y及びnormal shearで消せる。
E0 nonzeroではbase shearsの係数matrixのdetは18E0^2 nonzero。
E0=0の補集合ではb=-8、A0=1で、base q-shearがC0を動かし、
normal shear Delta=(-3x+6y)partial_zがD0を動かすため両者を消せる。
このzero補集合を分母E0で落とさない。
全模型は上のE0、k20=E0(b+8)、

\[
 Z=(x+3E_0y)\partial_y-3E_0z\partial_z,
 \quad a\ne4,\quad a+b+4\ne0
\]

という全二parameter源となる。

## 8. 残る全sourceと実行判断

constant-root tangentは三fixed P0、triangleは別noteの二parameter/oneparameter/twofixed、
concurrentは本稿の四/三/二parameter源へ還元された。
moving rootはtangentのoneparameter P0だけに尽くされた。
transverse conic、triangle h=L^3、全moving split qの全空性を別に保持する。

これらは必要全leading modelsであり、allhigher kernelの粗和を正しいG boundと混同しない。
次のfinite jobを行うなら、wholeparameter sourceと単位零補集合を最初から固定し、
actual normalizer及び全30/9 Jacobiを使う最弱命題へ限定する。
nonreduced R、m3及び他surfaceは依然として別義務である。
本稿の新手算術/normalizerは独立読解前で、未受理の証明として保存する。
