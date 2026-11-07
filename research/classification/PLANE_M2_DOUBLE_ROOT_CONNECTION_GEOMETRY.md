# rank二plane・m2 double normal rootの全接続幾何

2026-11-07。主義務G/E。実Dは既約被約quintic、S=P2はrestriction rank二、
zero cubic R=(f=0)は非零reducedとする。S沿い位数は二、最初のnormal quadraticは
generic double root、すなわちG2≠0かつdiscriminant零と仮定する。
本稿は全leading/sourceの必要幾何を証明する。全higher sourceの次元又は帰属は未完了。
新しいJacobi profile、minor、GB又はjetの計算は行わない。

## 1. 全primitive factorizationは二degreeだけ

Hamiltonian associated gradedによりgamma=G2/fはnormal flat connectionのparallel tensorである。
Gauss/UFDと二gradingから

\[
 G_2=h\ell^2,\qquad \ell=a u+bv,\quad\gcd(a,b)=1
\]

と書け、a,bは同じbase degree d、hのdegreeは3-2dである。
従ってd=0又は1だけである。primitive rowの共通零点は保持する。
非primitive境界では共通因子をhへ戻す。hのdivisorはRに支えられる。
実際R外でconnectionはregular、parallel nonzero tensorはそこで零になれないからである。
d=1ならa,bはindependent linear formsで、p=(a=b=0)はR上にある。

conormal covector connectionをbetaと書く。beta(E)=I、f betaの各one-form係数はdegree二、
各entryはR沿い高々単極であり、flatness conventionはd beta-beta wedge beta=0とする。
primitive root lineは保存され、

\[
 d\ell+\ell\beta=\theta\ell,\qquad d\theta=0.
\]

primitive rowのlocal unit成分で解けばthetaはR沿い高々単極。
closednessからthetaはconstant-residue logarithmic formである。
同様にt=tr betaはclosed logarithmic。parallel gammaのdeterminant twistを含めると

\[
 2\theta-t=d\log(f/h).\tag{1}
\]

Euler theta=d+1、Euler t=2であり、式(1)のdegreeとも一致する。
log性を全matrix betaへ課してはいない。

## 2. moving linear root：固定traceでprojective connectionは一意

d=1、q=f/h、ell=(a,b)とする。同じtとell,hを持つ二つのconnectionの差delta betaは
ell delta beta=0かつtr(delta beta)=0なので、有理one-form omegaにより

\[
 \delta\beta=\omega N_\ell,\qquad
 N_\ell=\begin{pmatrix}-ab&-b^2\\a^2&ab\end{pmatrix}.
\]

Nellのentriesのgcdは一。f delta betaはpolynomialかつ係数degree二であるため、
f omegaの各係数はpolynomial degree零、すなわちconstant。
一方delta beta(E)=0より(f omega)(E)=0であり、constant-coefficient one-formのEuler contractionが
恒等的に零になるのはzero formだけである。従ってdelta beta=0。
これは全global connectionの一意性であり、pointwise stabilizerが零という主張ではない。

## 3. moving rootのbasepointはqにも属する

pはR上にある。pがh上かつq外にあると仮定して矛盾を導く。
局所座標をh=x、ell=(x,y)とし、tracefree partを

\[
 B=A\,dx/x+C\,dy/x
\]

と書く。A,Cはregular matrices、trace零。regular部分もこの表示に吸収する。
delta=theta-t/2=(1/2)dlog qはp近傍でregularであり、root equationは
d ell+ell B=delta ellである。
x倍してx=0へ制限すると、A,Cのsecond rowsは零である。
trace零も使うと、両residue matricesはstrict upper triangular。
flatnessのx^-2係数はC|x0と[A|x0,C|x0]の符号付き等式であり、後者は零。
従ってC|x0=0、接方向dyのpoleは消える。
dyのroot equationをpで評価すれば、d ellのsecond成分は一、他はell(p)=0から零で矛盾する。
よって全sourceでpはq上にある。

逆にq(p)=0ならqはideal(a,b)に属し、linear c,dを選んで

\[
 M=\begin{pmatrix}a&b\\c&d\end{pmatrix},\quad\det M=q
\]

と出来る。second rowは補助frameなのでprimitive性を要求しない。
全completionの差はsecond rowにconstant倍first rowを加えるだけであり、M^-1dMは変わらない。
任意のallowed closed logarithmic trace tに対し、

\[
 \theta=(t+d\log q)/2,\qquad
 \beta=\theta I-M^{-1}dM\tag{2}
\]

はsimple-pole global source、Euler beta=I、flatであり、最初のrow equationを満たす。
前節の一意性から全moving-root connectionsが式(2)で尽くされる。
ただし、これはconnectionの存在であってPoisson P0又はfullhigher sourceの存在ではない。

## 4. actual P0によりbasepoint位置を更に制限する

全normal filtrationのweight零tensorP0はPoissonであり、その実五次式は正確にh Q^2、
Q=a u+bvはrank四quadric、vertexはpである。
従って既済の `DOUBLE_QUADRIC_SUPPORT.md` 第4節
「rank四Q、Hがvertexを避ける全配置は空」をP0へ適用出来る。
これによりpはhにも属する。式(2)のconnectionだけなら存在したp in q minus hをここで除く。

従って全moving-root Poisson sourceはp in h intersect qへ還元される。
reduced reducible cubicの四型ごとに有限marked配置である。
conic+lineではhは唯一のlineで、pはconicとの交点（transverse二点は同値、tangent一点）。
triangleではhは三線の一つ、pはhと他線の交点（全選択は同値）。
concurrent三線ではpは共通交点である。
ここでnormal右GL2及びmarked hを保つbase変換だけを使い、root-row mixingは使わない。
irreducible reduced cubic fにはline因子hが無いのでmoving-rootは存在しない。

Qのvertexを通るhはspan(a,b)の一線であり、Q restricted hは二linear factorsの積となる。
従ってこのP0は既済rank四H+2Qのcontact配置に入る。
そのseedpointのMPS帰属から未知の全P1/P2まで被覆されたとは推論しない。

## 5. constant root：全nonlog entryは閉有理one-formへ

d=0では定数normal GL2でell=uとし、

\[
 \beta=\begin{pmatrix}\theta&0\\b&\sigma\end{pmatrix},\qquad
 \theta-\sigma=d\log(f/h),\quad
 d\theta=d\sigma=0.
\]

flatnessのoff-diagonal条件はdb+(theta-sigma) wedge b=0なので

\[
 \omega=(f/h)b\quad\hbox{はclosed},\qquad
 \omega=N_2/h,\quad N_2=f b\in H^0(\Omega^1_{\mathbb P^2}(3)).\tag{3}
\]

これは全nonlog bを保持した正確な必要十分connection記述である。
thetaはRの全closed logarithmic periodsを保持し、sigmaは上の差で定まる。
Euler theta=Euler sigma=1、omega(E)=0。
hのdegree三とsupport条件から、conic+lineではh=f又はL^3、
three-linesではhはexponents111、21、3の全markedcomponent choicesである。
Rごとの同値を使う前にhの全選択を保持する。
式(3)のclosed rational formの全解を更にnormal-form化すること、
actual P0の存在とnormal bracket、全higher Jacobiを解くことは次の義務である。

## 6. 既約actual Dに必要な一つの全higher係数

constant rootのcanonical tracefree mixed gaugeではP0のu-columnはu倍であり、
normal bracketのv^2項はtarget G2 restricted u0=0から零。
全P1でu-rowのuで割れない項はv^2 C（Cはconstant base vector）だけ、P2にはnormal columnが無い。
C=0なら全u-rowがu倍となりactual Dはu因子を持つので、指定locusではC nonzeroが必要である。

P0のdiagonal mixed fieldsをX,Y、normal uv coefficientをkappaとすると、
全projective Jacobiのnormalweight一係数及びbase curlから、squarefree性に依らず

\[
 (2\kappa I+4DY-DX)C=0\tag{4}
\]

が必要である。原導出は `PLANE_M2_NILPOTENT_DIAGONAL_ROOTS_EXCLUSION.md` 第2、3節。
同導出ではoff-diagonal field、他のP1及び全P2は当該係数へ寄与しない。
normal doubleの場合、uv target coefficientは零なので、squarefreeの場合のkappa nonzeroは使えない。
式(4)をscalar+nilpotent invertibleと呼んで全空性を結論しない。
ここまでを次のwhole-source geometryの入口とし、period resonanceとC nonzeroの全補集合を保持する。

## 7. 証明状態と残余

全factorization、moving-root connection一意性と存在、p in h intersect qへの必要還元、
constant-rootの全closed-form module及びbad-rowの必要条件は手証明である。
moving connection部分は別作者によるchatでの独立読解済み、正式receiptは別に保存する。
fullhigher sourceのE帰属又はG上界は未証明。
nonreduced R、normal order三以上、plane以外のsupportは全て別残余として保持する。
