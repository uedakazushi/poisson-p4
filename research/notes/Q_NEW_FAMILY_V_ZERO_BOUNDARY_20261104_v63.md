# v=0、u≠0のfull tensor：Ga stabilizer、24次元境界、incidenceの区別

2026-11-04。対象は元の実際の族P(u,1,v)のv=0部分で、u≠0を保つ。
uvp開集合の結果をN²=0へ外挿せず、新しい境界族の全係数を扱う。
N²=0だけから既知N成分への所属を推論しない。

## 1. 実際のdivisorとmodular型

元full Pから直接divergenceを計算すると、homogeneous modular matrixは

\[
N=\begin{pmatrix}
0&0&0&0&0\\0&0&0&0&0\\4u/3&0&0&0&0\\
1/48&u&0&0&0\\0&0&0&0&0
\end{pmatrix}.
\]

N²=0、rank N=2で、ImN=span(ey,ez)のunit minorは4u²/3。
actual anticanonical divisorは引き続き3H+Qで、保存した5×5 quadric matrixのdetは

\[
\det Q=32768u^4/531441\ne0.
\]

従って全u≠0でQは滑らかで幾何的に既約。Hは唯一のtriple hyperplaneで、Q|Hのvertexは
同じz方向である。generic uvpのN²≠0との相違は実際のfull affine invariantであり、
この族を同じuvp orbitへ入れることはできない。

## 2. 全u≠0のstabilizerを尽くす

H、vertex、ImNから全affine旗を得る。ここでもPoisson scalar κを未知とする。
一般の低bracket比較とNyからη=κf=1、b=64c、d=a、D=a²−64c²、f=D、κ=1/D。
さらにN M=κ M Nのz,w成分はu b=0となるため、u≠0からb=c=0。
z,x成分はj=a³、従ってf=a²、κ=a⁻²となる。aは単元である。

この段階のy変換はex+a²y+gw+Tで、e,gを先に零としない。
full Pyz二次係数は

\[
[w^2]:a^3-a^2-64ag,\qquad
[xw]:-a(128e+ug).
\]

二次条件だけではu²=256に見かけの例外を作り得る。
しかし**full Pzwの一次w係数は正確に−2g/3**である。
よってg=0、次いでw²条件からa=1、xw条件からe=0。
これで全u≠0を尽くし、共鳴値を除外する必要がない。

Pxy、Pywのconstant coefficientとNzのconstant coefficientはtx=tw=i=0を与える。
残るz変換をz+Hx+Kw+Z、y平行移動をTとすると、三つのconstant条件は

\[
\frac89K=0,\qquad \frac89H=\frac23T,\qquad
-H/36+4uZ/3=0.
\]

係数matrixのdetは−256u/243というunitで、唯一の解は

\[
H=3T/4,\quad K=0,\quad Z=T/(64u).
\]

これを元full六成分の全degree0,1,2係数へ戻して恒等零を確認する。
全homogeneous行列はI+3TN/(4u)そのもので、N²=0だからこれはexp(3TN/(4u))。
従って

\[
\boxed{\operatorname{Stab}_{PGL_5}([P(u,1,0)])=\mathbb G_a
\quad\text{for every }u\ne0.}
\]

有限E成分が境界にも残ると仮定してはいけない。uvpでのE表示は1/vを持ち、
上の全係数比較がこの境界で追加有限成分の不存在を示す。

## 3. parameter orbitと飽和の次元

source(u,0)、target(U,0)についても全affine旗とNyが有効で、低bracketはη²=1、U=ηuを
強制する。旧full tauがu→−uを実現するので、全u,U≠0で

\[
P(u,1,0)\sim_{PGL}P(U,1,0)\iff U^2=u^2.
\]

parameter quotientはGm上のu→−uによるGm_A、A=u²である。
各orbitはdim PGL5−dim Ga=23。parameter同定は有限なので、この一変数族の
PGL saturationの閉包は既約**24次元**となる。
元の三parameter族のv→0 specializationなので、既受理の新25次元成分の境界に含まれる。
それ自体を新しいirreducible componentと呼ばない。

## 4. 実際の旧99-variable sourceへの所属

ここはN²=0だけによる推測ではない。凍結された
Q_POSITIVE_WEIGHT6_HQUARTIC_EXCLUSIONのexact99 polynomial graphへ
(n0,n1,n2)=(u,1,0)を代入した全99値を保存する。
最初の90係数を元positive Pへ戻すと、今回のP(u,1,0)の全六式と厳密に一致する。
元209式への全pullbackを同じ代入で確認すると全て零になる。
唯一の非零一般pullbackだったrow179=(4/3)n0²n2も零である。

したがってこの族は**実際のclosed H/compact-center/projective-pencil/N(P)² incidence**に入る。
しかしそのmarkは

* center=[0:0:0:1:0]で、有限centerではない。
* constant pencilは∂y∧∂zで、非退化ではない。

全九個のauxiliary chart座標が零である。このclosed incidenceへの所属を、
finite-center・非退化pencilの部分の閉包への所属と同一視しない。

## 5. 二種類の「N」と既知成分を混同しない

旧99 sourceのNはaffine modular matrix N(P)である。
CAとB18で使う既知27次元の**MPS N normalization family**とは別物である。
MPS Nのproper固定旗19空間はCA20空間内のf∈hWという部分空間であり、N⊂CA。

独立受理済みのglobal proper CA/CB模型のactual divisorはそれぞれ
g h³ f(h,x)、h·a b F(h,x)²である。幾何的には全てlinear factorだけの五次式、または零。
proper模型の像が成分全体なので、閉包の特殊点もこの条件を保つ。
今回のnonzero divisor3H+smooth irreducible Qはその型でない。従って全u≠0で

\[
P(u,1,0)\notin CA\cup CB\cup N_{\rm MPS}.
\]

旧R4_H_QUARTIC_ESSENTIAL5も、この族の所属先とは結論できない。
むしろその既受理componentはH+irreducible quarticを持つ一つの23次元orbitの閉包である。
今回のorbitも23次元だがdivisor型が違うのでそのopen orbitには入らず、orbit閉包の境界は
次元≤22だからboundaryにも入らない。従ってこのHquartic componentからも全u≠0を排除できる。
closed209-incidenceがそれより大きなboundary成分を持つことと矛盾しない。

## 保存と範囲

専用code `q_new_family_v_zero_boundary_20261104_v63.py` と同名report directoryへ
全P、N、Q、決定的な一次/二次係数、unit系、Ga行列、全99入力と209pullback確認を保存する。
新しい点・素数・GB・非線型同値solverは使わない。

```sh
source /workspace/cloud-setup/activate.sh
.venv/bin/python code/q_new_family_v_zero_boundary_20261104_v63.py \
  --output work/vzero-boundary-replay
```

対象はこの明示的なv=0族とその飽和だけで、u=0、p=0、whole新成分の全boundary、
その他のknown componentへの所属・非所属は分類しない。専用独立監査前は採用を保留する。
