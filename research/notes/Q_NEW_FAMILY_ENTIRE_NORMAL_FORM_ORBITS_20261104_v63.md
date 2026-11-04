# uvp normal form全域のPGL同定と有限parameter quotient

2026-11-04。対象は実際の全Poisson normal form P(u,1,v)の領域

\[
\mathcal U=\{u v(32u+v)\ne0\}
\]

全体である。source、targetの両方をこの領域に置く。genericなlogarithmic residueや
invariant lineの一意性を使わず、低次のfull bracket係数だけから必要条件を出す。
これによりu=±8、u=±16、J=0を含む全parameter同定を尽くす。
normal formのPGL saturationのorbitを分類するが、全C_NEW25のsmooth-Q点がこのnormal formに
属するとは主張しない。

## 1. 全域で正当なaffine旗

実際のanticanonical divisorは3H+Qで、uvp領域ではQは滑らかでHは唯一のtriple
hyperplaneである。Q|Hの唯一のvertexはz方向。さらにdivergenceから得るhomogeneous Nの
像は、minor4u²/3により全u≠0でspan(ey,ez)となる。target(U,V)でも同じである。

従って任意のfull PGL同値はHを保つaffine変換であり、vertexとImNも保つ。
Poisson scalar κを未知のまま残すと、その全一般形は

\[
x'=ax+bw+t_x,\quad w'=cx+dw+t_w,\qquad
y'=ex+fy+gw+t_y,\quad z'=hx+iy+jz+kw+t_z,
\]

で、f j(ad−bc)≠0、κ≠0を満たす。u=±8でcompatible constant bivectorの空間が跳ぶことは
既受理だが、ここではその一意性に依存して旗を作らない。

## 2. sourceとtargetのuが±でしか違わないこと

恒等式κ A P(u,1,v) Aᵀ=P(U,1,V)(Ax+t)を比較する。
full bracketの三成分

\[
P_{xw}=8/9,\qquad P_{xy}=1/36-256w/3,\qquad P_{yw}=4x/3
\]

はu,vによらない。constant Pxwからκ(ad−bc)=1を得る。
eta=κfと置く。Pxyのx,w係数とPywのx,w係数から

\[
\eta b=64c,\quad\eta a=d,\quad
\eta d=a,\quad\eta c=b/64
\]

となる。aやcを先に割ることはしない。実際

\[
\eta(\eta a-d)+(\eta d-a)=(\eta^2-1)a,
\]
\[
\eta(\eta b-64c)+64(\eta c-b/64)=(\eta^2-1)b.
\]

x/w blockは可逆なのでa,bは同時に零ではない。従ってeta²=1である。

一方、affine divergence covarianceは

\[
N_{\rm target}M=\kappa M N_{\rm source}
\]

である。y,h成分は4U/3=κf·4u/3を与えるので

\[
\boxed{U=\eta u=\pm u}.
\]

ここまでu,v,U,Vが表示されたdomainにあること以外のgeneric条件を使っていない。
既知のfull tau:(U,V)→(−U,−V)で符号を変えればU=uへ戻せる。

## 3. 二本の直線の変換則と全域でのJ

H上ではtranslationが消える。l±=x±8wと置くと、上の低bracket条件だけから次を得る。

* eta=1ではb=64c,d=a、l±targetのpullbackはそれぞれ(a±8c)l±source。
* eta=−1ではb=−64c,d=−a、l+targetのpullbackは(a+8c)l−source、
  l−targetのpullbackは(a−8c)l+source。

det(x/w block)=eta(a²−64c²)は非零なので、二つのscale(a±8c)は共に非零。
従って**与えられたnormal formの二本のlineは、全ての実際の同値でunordered pairとして保たれる**。
特殊parameterでfoliationに他のinvariant lineが現れる可能性を排除する必要はない。
genericなline一意性の定理を全域へ外挿していない。

vertexからのP2 quotientで、Q|Hのconicはscalar倍して

\[
q_2(u,v)=-x^2+64w^2+64uxy+32uvy^2
\]

である。実際のfull projective同値はHとQを保つので、このprojected conicも非零scalar倍まで
保つ。上記line covarianceと合わせ、dual-conic pairingの比

\[
J=256\frac{(\ell_+^{\mathsf T}q_2^{-1}\ell_-)^2}
 {(\ell_+^{\mathsf T}q_2^{-1}\ell_+)(\ell_-^{\mathsf T}q_2^{-1}\ell_-)}
\]

も保たれる。q2はその対称行列を同じ記号で表した。
これはmatrix inverseのprojective covariance、各lineのrescaleと交換に不変な比である。
既受理の全有理式は

\[
\det q_2=-2048u(32u+v),\quad
\ell_\pm^{\mathsf T}q_2^{-1}\ell_\pm=32u/(32u+v),
\qquad J=(v/u+16)^2.
\]

全分母とself pairingはuvp領域全体で非零である。従ってsourceとtargetについて

\[
\boxed{(V/U+16)^2=(v/u+16)^2}
\]

が必要である。これもu=±8、u=±16、J=0を含む。

## 4. 全parameter同定：必要十分性

U=±uと二乗の等式から、targetの可能性は正確に

\[
(u,v),\quad(-u,-v),\quad(u,-32u-v),\quad(-u,32u+v)
\]

の四つに限られる。旧full tauと、独立監査済みのfull sigma行列は表示されたuvp領域全体で
定義され、その四つの候補を全て実現する。よって

\[
\boxed{P(u,1,v)\sim_{\rm PGL}P(U,1,V)
\iff U^2=u^2,\ (V/U+16)^2=(v/u+16)^2}
\]

が**normal form領域全体**で成立する。二つのnormal formの係数が同じ標数0のfield kに
属していれば、同じ式から符号もk内の±1であり、実現matrixもk上有理な式なのでPGL(k)で
同値である。一般のmoduli k点が必ずk-valued normal formを持つという別の主張はしない。

J≠0では四つのnormal formは相異なる。J=0ではv=−16uでsigmaがparameterを固定し、
u≠0なので相異なるnormal formは二つだけとなる。追加parameter同定が特殊uで現れる余地はない。
同一tensorのstabilizerの分類とは別の命題であり、そちらは専用whole-open packetを使う。

## 5. finite parameter quotientのcoordinate ring

z=v/u+16と置く。uは単元なので、parameter ringは厳密に

\[
k[u,u^{-1},z,(z^2-256)^{-1}]
\]

である。tauはu→−u、z固定、sigmaはu固定、z→−z。従って二つの独立なsign作用の
不変部分は、分母が既に不変であることとmonomialの偶奇から

\[
\boxed{k[A,A^{-1},J,(J-256)^{-1}]},\qquad A=u^2,\quad J=z^2.
\]

元ringはこの不変ring上でbasis1,u,z,uzのfree rank4である。
従って有限parameter quotientは

\[
\boxed{\mathbb G_{m,A}\times(\mathbb A^1_J\setminus\{256\})}.
\]

幾何的にJ=0のfiberは二つのu-rootとz²=0からなるlength4であり、四つの被約点ではない。
orbit内の相異なるnormal form数が二つという主張と整合する。
J=256はv=0または32u+v=0に正確に対応するので、今回のdomainから除かれる。

このfinite-group quotientはparameter空間のaffine quotientとして実際に計算したもので、
whole componentへのgood/categorical PGL quotientのregularityを仮定していない。
このquotientの幾何点は、displayed normal formのPGL saturationに属するorbitを分類する。
そのsaturationがwhole component内の全smooth-Q点を尽くす、あるいはこのparameter quotientが
whole componentのglobal moduli schemeである、とは主張しない。

## 保存と確認

専用producerは旧gateとhash、四sign作用の有限な組だけを確認する。新しいCAS、点、素数、
GB、不変式評価はない。source codeは
`code/q_new_family_entire_normal_form_orbits_20261104_v63.py`、同名report directoryに
全低係数関係・line covariance・有限quotientのproof schemaを保持する。

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_entire_normal_form_orbits_20261104_v63.py \
  --output work/entire-normal-form-orbits-replay
```

rootの独立概念監査へ渡してから採用する。uvp=0の境界とwhole componentの未認識点は対象外。
