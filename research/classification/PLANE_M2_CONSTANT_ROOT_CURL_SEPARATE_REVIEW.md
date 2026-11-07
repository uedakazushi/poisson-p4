# m2 constant root：全SNC曲線とnilpotent-curl空性の独立検証

2026-11-07、独立 reviewer `/root/quartic_review`。**PASS（純証明）**。
新しいsource/minor/GB/jet/数値Jacobi計算はない。
`PLANE_M2_CONSTANT_ROOT_MPS_AND_NILPOTENT.md`（d1dba303）の**§§1–3だけ**、
`PLANE_M2_NILPOTENT_DIAGONAL_ROOTS_EXCLUSION.md`（4a27cbcd）の全raw係数証明を独立導出した。
前者§§4–5のflag02専用rational-affine別証明は今回の合成に不要で、独立受理の依存にしない。

対象は rank2 restriction plane、reduced非零zero cubic R、actualD被約既約quintic、
normal位数m2、generic squarefree normal quadratic、primitive degree flags(0,q)、q=0,1,2の**全源**。
normal double-root、R非被約、m3等は対象外。
定数normal GL2及びprimitive factorizationから

```
G2=κ0 h u(Aq u+Bq v),  h Bq=f,  κ0≠0,
M=[[1,0],[Aq,Bq]],  gcd(Aq,Bq)=1
```

となる。Aq全係数、全logarithmic residues、全higher sourceを保持する。
row共通因子境界はhへ戻し別flagとして保持、q3は既root-line class矛盾により空。

## 完全なnormal-weight零源

constant root uのparallel性からβuv=0。
flatnessによりdiagonal entriesはclosed、R被約のsimple-pole boundからconstant-residue logarithmic。
off-diagonal βvuのnonlog項は除かない。
base Jacobian代表、tracefree mixed Euler gaugeを取ると

```
P0_ij=(P_f)_ij,  P0_i,u=u X_i,
P0_i,v=u Z_i+v Y_i,
P0_u,v=k20 u²+κ uv+k02 v².
```

X,Y,Zは全linear base fields。target G2|u0=0とactualEuler-Pfにより3f k02=0、k02=0。
このgaugeは全normal-linear Euler補正をK内へ吸収する正当な全源表示。
完全P1はraw36/Euler6の30、P2は全9で、固定切片ではない。

## triangle及びtransverse conicの全actual SNC曲線

triangle f=xyzではbase weights0、(wu,wv)=(N+1,N)、N≥2。
P0のbase/diagonal mixedはweight0、off-diagonal u→vとK20u²はpositive、Kuvは0。
K02v²と逆mixedは無い。全raw36のP1もweight≥N−1>0、全P2は≥2N。
Aqやresidueのzero/special値も同じ不等式に含む。

transverse conic f=z(xy−z²)ではdiagonal mixed fieldsは
`x∂x−y∂y` のconstant倍。これは全closedlog diagonal formの
`dlog(xy−z²)−2dlogz` に対応する全horizontal spaceである。
weights(base;u,v)=(0,0,1;N+2,N)、N≥4を取る。
baseJacobian各項は非負、diagonal mixedは0、arbitrary linear Zのoff-diagonal項は≥1、
K20は2、Kuvは0。全raw36 P1は≥N−3>0、全9 P2は≥2N−1>0。
fullAutやnonlog Zを消してこの評価を得たわけではない。

従って両型で
`Pt_ij=t^(-wi-wj)Pij(t^w X)` はtの多項式、t≠0で元とPGL同値、t0でもprojective Poisson。
actualDの変換は `Dt=t^(-Σw)D(t^w X)`。
triangleではΣw=2N+1、G2のfu vだけがrelative0、hAq u²はpositive。
transverse conicではΣw=2N+3、fu vの唯一最低項はxyzuv、hAq u²は≥1。
normalorder r≥3項はそれぞれ≥(r−2)N−1、≥(r−2)N−3で全て正。
**actualD limitはκ0xyzuv、κ0≠0**で真のSNC log-symplectic点。

このactual irreducible曲線全体を含むambient componentはSNC locusに会う。
既受理MPS Theorem1.4の**SNCに会う全component**scopeでそのcomponentはMPS40、元の任意源も入る。
一点のknown closed-union membershipを曲線へ逆向きに外挿したり、T01一意性を用いない。

## tangent/concurrentの全nilpotent diagonal space

tangent conic f=z(xz−y²)ではX,Yは `2y∂x+z∂y` のconstant倍。
concurrent f=xy(x−y)ではX,Yは全 `(a x+b y)∂z` space。
全residue/zero valuesを含め、traceX=traceY=0、Xf=Yf=0、X∧Y=0で、
4DY−DXの全linear combinationはnilpotentである。
actualG2のuv係数は3fκ。X∧Y/Xf/Yfのcross termsは零、Zはu²項だけ。
targetはκ0f≠0なので**κ≠0**。κ0=0はm2を失う別義務として保持する。

## projective Jacobiから一つのbad coefficientを全域排除

全P1のu-columnのu非可除部分は `v²C`、Cはconstant base vector。
P2はnormal columnsを持たない。他の全higher coefficientsは固定しない。
strict Jacobi自体を零と仮定せず、実際のprojective条件 `E∧Jac(P)=0` を使う。

u=0でnormal v³/v²係数を直接並べると

```
Jac_iuv[v³]=κ C_i+(DY C)_i,
Jac_iju[v²]=[(2Y−X)∧C]_ij−(P_(C·∇f))_ij.
```

一式目はP_i,u ∂uP_uvとP_u,l ∂lP_v,i、
二式目はmixed-vの2v derivative、mixed-uのu derivative、及び−C·∇P_fだけである。
Z/k20はu因子を持ち、他P1/P2の項は別normal monomial/高weightなので寄与しない。
projectivequad(i,j,u,v)のv³係数は**正確に**

```
Ebase∧(κC+DY C)−(2Y−X)∧C+P_(C·∇f)=0.
```

base volumeでbivectorをvectorへ双対化しcurlを取る。
`curl(E×D)=−2D`（D constant）、
`curl(V×C)=DV C`（V tracefree linear、C constant）、
`curl(grad(C·∇f))=0`。
従ってglobal signを整理して

```
(2κ Id+4DY−DX) C=0.
```

同じ式はcolumn divergenceを使っても得る。sourceのslice又はstrictJacobiへのすり替えはない。
scalar2κ≠0とnilpotent matrixの和は可逆なのでC=0。
X=Y=0等のspecial periodも同じ議論に含む。
全u-rowがu倍、normal bracketもu(k20u+κv)、P2にnormal行は無い。
従ってu0ではprojective bivectorがhyperplaneの三次元tangent spaceに入りrank≤2、
actualDがuで割れる。actualD既約に反し、両nilpotent R型の全constant-root源は空。

## 全reducedR m2 squarefreeへの合成

新constant-root結論はtriangle/transverse-conicで全MPS、tangent/concurrentで全空。
既受理root-line lemmaでflags03/12は空、conic/splitのflag11は全MPS又は空。
既約被約Rの全rank2plane排除も既受理の全smooth/node/cusp theoremをそのまま適用する。
これで全reduced Rと全六degree flagsを尽くし、指定**全m2 squarefree sourceはMPS40閉和に入る**。

aggregate `PLANE_RANK_TWO_M2_SQUAREFREE_REDUCED_CUBIC_COVERAGE.md`（920d2580）の全表を照合した。
split旧metadataはVM reset後exact復元できず、数学noteは旧hash同一のまま
`...SPLIT_SEPARATE_REVIEW_RECOVERY_V2_20261107.json` を依存とする。
conic/root-line receiptsは旧bytes同一。新計数は無く、全E/G又は他支持の完了は主張しない。
normal double quadratic、R非被約、m3、他surfacesは残す。
