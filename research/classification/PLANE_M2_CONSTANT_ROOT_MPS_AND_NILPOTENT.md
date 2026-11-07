# m2 constant normal root：全SNC認識とflag02のnilpotent排除

2026-11-07。主義務G/E。全sourceの手証明、計算jobなし。
実D被約既約quintic、rank二restriction plane S、reduced nonzero zero-cubic R、
S沿い位数二、normal leading quadratic generic squarefreeとする。
root-matrix degree flags(0,q)の全源を扱う。q=0,1,2で、q=3は既line-class矛盾で空。

**結論。** Rがnonconcurrent triangle又はtransverse conic+lineなら全constant-root sourceはMPS40。
flag(0,2)でRがtangent conic+line又はconcurrent triangleなら全sourceは空。
従ってreducedRのflag(0,2)は全てMPS又は空となる。
tangent/concurrentのflags00/01、normal double quadratic、nonredR及びm3等は明示残余である。

## 1. 全triangular connectionとraw polynomial source

primitive root factorizationを定数normal GL2で

\[
 G_2=\kappa h\,u(A_q u+B_q v),\quad hB_q=f,\quad\kappa\ne0,
 \quad M=\begin{pmatrix}1&0\\A_q&B_q\end{pmatrix}
\]

と取る。Aq,Bqは同じdegree q、gcd(Aq,Bq)=1を保つ。
q=0はconstant root frameでA0=0,B0=1と出来る。
全line-class証明によりθ1,θ2はR沿いconstant-residue closed logarithmic formsで、

\[
 \beta=M^{-1}\operatorname{diag}(\theta_1,\theta_2)M-M^{-1}dM
 =\begin{pmatrix}
 \theta_1&0\\
 [A_q(\theta_2-\theta_1)-dA_q]/B_q&\theta_2-d\log B_q
 \end{pmatrix}.
\]

diagonal entriesはclosed logarithmic、off-diagonalはnonlogを含む全源である。
βuv=0なのでcanonical tracefree mixed gaugeでP0のu-columnはu倍、
v-columnの唯一のoff-diagonal項はu倍である。
target G2|u0=0からactualPfaffian coefficient3f k02=0となり、raw normal bracketにはv²項が無い。
K=k20u²+k11uvの全二係数はここで固定しない。

全higher sourceはP=P0+P1+P2、既受理の全30/9方向を保持する。
以下のweightsはAqの全coefficients、全logarithmic residuesと全P1/P2に一様であり、
first symbolのchosen pointだけを扱うものではない。

## 2. nonconcurrent triangleの全constant-root源はMPS

f=xyzとする。closed diagonal logarithmic formsから得るtracefree linear base fieldsは全てdiagonal。
base weightsを零、normal weightsをwu=N+1,wv=N、N≥2とする。

P0のbaseとdiagonal mixedはweight零、off-diagonal u→vはweight一、
normal Kのu²はweight一、uvは零で、反対向きやv²項は無い。
P0は全weight非負。
全三十P1方向の最小weightはN−1、全九P2方向は2Nで、全て正。

GL5から得るactual curve Pt=gt−1 P(gtX)gt−Tはtのpolynomialであり、
t≠0では元PとPGL同値、t=0もPoisson。
sum(weights)=2N+1に対し、G2のAu²部分はpositive relative weight、
fuvだけがweight零。normal order三以上のD項もrelative weight≥N−1>0。
ゆえに

\[
 D(P_{t=0})=\kappa xyzuv
\]

で真のSNC log-symplectic point。
irreducible curveを含むambient成分はSNC locusと交わるため、
既受理のMPS Theorem1.4の全component scopeによりMPS40の一つとなる。
元の全Pもその閉包に入る。
T01の一意性やSNC一点の閉和所属だけを使った推論ではない。

## 3. transverse conic+lineの全constant-root源もMPS

f=z(xy−z²)とする。
closed diagonal logarithmic formsのhorizontal部分はdlog(xy−z²)−2dlogzのconstant倍であり、
対応するtracefree base fieldsはx∂x−y∂yのconstant倍。
ここでfullAut(conic+line)やnonlog off-diagonal parameterを除いていない。

weights(x,y,z,u,v)=(0,0,1,N+2,N)、N≥4を取る。
base Jacobian tensorの各項はweight非負、diagonal mixedは零。
off-diagonal u→vは少なくともweight一、Kのu²は二、uvは零である。
全P1最小weightはN−3>0、P2は少なくとも2N−1>0。
従って同じactual GL5 curveがt=0までregularに延びる。

sum(weights)=2N+3に対して、fuvの最小項はxyzuv。
hAq u²部分のrelative weightは少なくとも一、normal order三以上も少なくともN−3>0。
従ってlimit actualDはκxyzuv、全sourceが同じMPS全component theoremへ入る。
q=0,1,2の全primitive rootsと全zero coefficient boundariesを保持した証明である。

## 4. flag(0,2)をactual h=∞ chartでscalar化する

q=2ではhは一つのlineであり、actualD(P0)=κh u(A2u+B2v)のsimple componentである。
line-class relationはθ1−θ2=−2dlogh。
rational normal change

\[
 u'=u,\qquad v'=(A_2u+B_2v)/h^2
\]

で両root-line connectionはθ1となる。
これはwholePのPGL同値とは主張せず、dense open上の必要恒等式を読むためだけのrational座標である。

h=1 affine chartを取る。hがP0のPoisson-invariant hyperplaneなのでP0のaffine coefficientsはquadratic以下。
一般homogeneous quadratic representativeから生じるcubic部分は、h-rowmodhがEuler倍であることにより相殺する。
したがってnormal-weight零bracket{u,v}のu² coefficientはbase-independent constantである。
ここで任意Euler代表のdivergenceやstrict stabilizerを仮定しない。

affine base PoissonはπS=3B∂r∧∂s、A=A2|h1、B=B2|h1。
scalar connectionにより{ri,u'}=u'Xi、{ri,v'}=v'Xiで、Xはθ1のHamiltonian contraction。
両mixed vectorsが同じなので、actual Pfaffianから

\[
 \{u',v'\}=k u'v',\qquad k\ne0
\]

を得る。k非零はD2非零、すなわちm二条件そのものであり、零補集合はm三以上の別義務。
以下の二nilpotent場合ではX(B)=0。
v'=Au+Bvをbracketに代入すると

\[
 \{u,v\}=\frac{kA+X(A)}B u^2+k uv.
\]

原affine bracketがquadraticなので、あるconstant cに対し

\[
 X(A)+kA=cB. \tag{1}
\]

が必要である。

## 5. tangent/concurrentのflag02はnilpotencyで全空

tangent conic+lineではh=z、B=x−y²。
θ1|h1はdlogBのconstant倍なので、Xは∂y+2y∂xのconstant倍である。
これはdegree≤2のpolynomialsの有限空間でnilpotent、Bを保存する。

concurrent triangleではhをその一成分に選び、affine base coordinatesをr,sとしてB=r(r−1)。
θ1|h1はdlogr,dlog(r−1)のconstant linear combinationなので、
X=(ar+b)∂s、やはりdegree≤2でnilpotentかつX(B)=0。

V=polynomials_degree≤2 / C B上でもXはnilpotentである。
式(1)は(X+kI)[A]=0、k≠0なのでX+kIはinvertible。
従ってA∈C Bとなる。
これはhomogenization後もA2=cB2であり、primitive row(A2,B2)に反する。
その非primitive境界はfactorization時にflag00等へ移すので、無断で別sourceを捨てていない。
よって両R型のflag02全sourceは空。
nilpotent X=0のspecial residuesも同じinvertibilityに含まれる。

## 6. 完了範囲と独立確認

全reducedRのflag02はMPS又は空へ尽くされた。
triangle/transverse conicのflags00/01も全MPS。
flag11は別二proofでMPS又は空、flags03/12はline-class矛盾で全空。
tangent/concurrentのflags00/01、generic double quadratic、nonredRとm3は未完了。
normal coefficientsのgeneric sample、minor、Jacobi profile、GB又はjetは一切計算していない。
本noteの正式独立receiptは別作者の手読解後に保存する。
