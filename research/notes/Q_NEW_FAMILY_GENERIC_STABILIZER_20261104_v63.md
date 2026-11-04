# NEW族のgeneric full PGL stabilizer：Ga × μ2

2026-11-04。元の実際のP(u,1,v)の**geometric generic projective stabilizerを尽くす**。
既知のmodular Gaに加え、二つのsigma変換の差から非自明なinvolutionを構成する。
さらにintrinsicな旗を保つ一般affine変換の全係数を比較し、残るgeneric成分がないことを示す。
特殊parameterでのstabilizerの跳びは分類しない。

## 1. 二つのsigmaからinvolutionを作る

受理済みのfull sigma行列をM+、M−とする。両方のPoisson scalar kappaは同じなので

\[
R=M_-^{-1}M_+
\]

は元のP(u,1,v)そのものを保つ。既知のhomogeneous modular行列をNとすると、全係数で

\[
RN=NR,\qquad R^2=\exp(t_0N),\qquad t_0=-1/(128uv)
\]

となる。N³=0なのでexp(tN)=I+tN+t²N²/2という有限多項式である。
従ってE=exp(−t0N/2)Rは、Nと可換なinvolutionである。簡約した全行列は

\[
E=\begin{pmatrix}
1&0&0&0&0\\
-4/(3v)&-1&0&0&0\\
4/(3v^2)&2/v&1&0&0\\
-2/(3v^3)&-1/v^2&0&1&0\\
1/1536&0&0&0&-1
\end{pmatrix}
\]

である。行・列の順は(h,x,y,z,w)。E²=I、EN=NEを直接検査する。
modular Gaのx,wのlinear blockはidentityだがEは−identityなので、EはGaに入らない。
これだけでは全stabilizerを尽くしたことにならず、次の係数比較を必要とする。

## 2. 全generic PGL変換を旗付きaffine形へ落とす

実際のanticanonical divisor3H+Qの唯一のtriple hyperplane Hを保つため、generic stabilizerは
Hの補集合でaffine変換である。さらにQ|Hの唯一のvertexはz方向であり、既受理の
unique compatible constant bivector lineはC=∂y∧∂zである。affine変換はconstant bivectorと
Schouten compatibilityを保つので、このlineとそのdirection planeも保つ。
従って、一般のstabilizerは次の形に尽くされる。

\[
x'=ax+bw+t_x,\quad w'=cx+dw+t_w,\quad
y'=ex+fy+gw+t_y,\quad z'=hx+iy+jz+kw+t_z.
\]

invertibilityはf j(ad−bc)≠0を要する。Poisson scalar kappaも未知として保持し、
kappa A P Aᵀ=P(Ax+t)を比較する。κ=1を先に仮定しない。

constant Pxwからκ(ad−bc)=1、PxyとPywの一次係数から、eta=κfについて

\[
\eta b=64c,\quad\eta a=d,\quad\eta d=a,\quad\eta c=b/64
\]

を得る。divergenceのaffine covarianceはN'=κ A_*Nであり、targetは同じPなので、
Ny=4u/3のconstant coefficientからeta u=uとなる。generic u≠0なのでeta=1。
これでb=64c、d=aとなる。D=a²−64c²と置くと

\[
f=D,\quad\kappa=1/D.
\]

二次Pxz,Pzwの共通因子x+vyを比較すれば

\[
j=D^2,\quad e=(D-a)/v,\quad g=-64c/v.
\]

係数因子を消す際にx/w blockの非退化性だけを使い、aまたはcが非零という一方のchartを
固定しない。Pxy,Pywのconstant coefficientとNzのconstant coefficientは

\[
t_x=-c/48+2(a-1)/(3v),\quad
t_w=(1-a)/3072+2c/(3v),
\]
\[
i=3D(t_x+vt_y)/4+D(1-D)/(64u)
\]

を与える。残る未知数はa,cとt_y、さらにh,k,t_zだけである。

## 3. 全quadratic caseを尽くす

p=32u+vとし、残るPyz二次係数を次で記す。

\[
Q_1=-pD^2+pa^2+4096ac-64vc^2,\quad
Q_2=-D^2+a^2+uac+64c^2,\quad
Q_3=vD^2-va^2+4096ac+64pc^2.
\]

実際のx²,xw,w²係数はQ1/(64v)、64Q2/v、Q3/vで、yを含む二次係数は恒等的に消える。
全必要条件はQ1=Q2=Q3=0である。次の全係数恒等式を保持する。

\[
Q_1+Q_3+(v-p)Q_2=-32ac(u^2-256),
\]
\[
Q_3+vQ_2=ac(4096+uv)+64c^2(p+v).
\]

generic u²−256≠0、p+v≠0では、最初の式からac=0、次の式からc=0となる。
D≠0なのでa≠0であり、Q2=a²(1−a²)から**a=+1または−1の二場合だけ**となる。
この議論は代数閉包上でも成り立つので、非有理なa,cの追加caseを残さない。
u²=256やv=−16uなどの特殊値へ同じexhaustionを外挿しない。

## 4. 残る全一次・定数係数

t_y=Tを自由に残す。どちらのaでも残る五つの係数条件はh,k,t_zの線形式であり、
certificateの3×3 minorが非零なので解は一意である。

* a=1では h=3T/4、k=0、t_z=3vT²/8+T/(64u)。全行列はexp(3T N/(4u))そのもの。
* a=−1では h=(3Tv²−8)/(4v²)、k=0、
  t_z=(72T²uv³−192Tuv+3Tv²−4)/(192uv²)。全行列は
  exp(3(T−4/(3v²))N/(4u))Eそのもの。

全六成分のdegree0,1,2係数への代入と、これらのmatrix等式を全有理式で検査する。
従って二つの発見されたmatrixだけをもって全群と呼ぶのではなく、intrinsicな旗から
始めた全generic係数caseをこの二familyで尽くした。

## 5. 群としての結論と範囲

EはNと可換で位数2、EはGaの外にある。従って

\[
\boxed{\operatorname{Stab}_{\operatorname{PGL}_5}([P(u,1,v)])
       \simeq\mathbb G_a\times\mu_2}
\]

がk(u,v)上およびその代数閉包上のgeneric stabilizerである。群同型は
(q,epsilon)→exp(qN)E^epsilonであり、各componentでy平行移動Tからqを戻せるので、
単なるgeometric pointsの個数の主張ではない。標数0ではこのalgebraic subgroupはsmooth。
必要条件は最終的にD=1、κ=1を強制するので、全generic projective stabilizerは選んだ
Poisson代表もscalar変更なしで保つ。

parameter間の同定を与える既知Klein fourと、同一tensorのstabilizer Ga×μ2は別の群である。
全special parameter、boundaryのgroup scheme、generic family normalizer全体の群構造は
今回の結論に含めない。

## 保存と再現

`code/q_new_family_generic_stabilizer_20261104_v63.py`と同名report directoryに、R、N、E、
全quadratic combination、両caseの全五線形式・rank3 minor・一意解・全matrixを保存する。

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_generic_stabilizer_20261104_v63.py \
  --output work/generic-stabilizer-replay
```

一回のbounded workflowは60秒/2GB内。新しい点・素数・GB・非線型equivalence solverを使わない。
generic case exhaustionとfull homogeneous actionの独立確認を受けてから結論を採用する。
