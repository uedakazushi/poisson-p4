# rank二plane・reduced R・m2 squarefree normal coneの全MPS被覆

2026-11-07。主義務G/E。全sourceの幾何学的尽くし。
実Dが被約既約quintic、S=P2⊂SingDがrestriction rank二のPoisson plane、
πSのzero cubic Rがreduced、DのS沿い位数が二、normal quadraticのgeneric二rootが相異なるとする。

**定理。** この全projective Poisson sourceは既知MPS40成分の閉和に入る。
表中の空性は全source、MPSは全pointwise認識である。
sourceが既知成分へ入ることを全像dimension≤24と言い換えない。
normal quadraticのdouble-root枝、nonreducedR、m3や他支持はこの定理に含めない。

## 1. 完全root factorizationと全六degree flags

Hamiltonian associated-graded補題によりγ2=G2/fはparallel。
twisted binary-quadratic discriminantはinvariantなのでΔ(G2)=c f²、c≠0。
UFDとGauss lemmaにより、全pointで

\[
 G_2=h\ell m,\quad\binom\ell m=M\binom uv,\quad h\det M=f
\]

と書け、二root rowsはprimitive homogeneous polynomial vectors。
そのbase degrees d1≤d2は

\[
 (0,0),(0,1),(0,2),(0,3),(1,1),(1,2)
\]

の全六flagに尽くされる。
root-row basepointsは残し、unused componentsも許す。
rowに共通因子が現れる境界はhへ戻すので、別flagとして全て含む。
f reducedだからhとdetMは互いに素なR-component productsである。
詳しい全source証明は `PLANE_M2_SPLIT_ROOT_MATRIX_ENTRANCE.md`。

## 2. 二つのflagsはline-class矛盾で空

primitive rowのunit成分からroot-line connection θiはR沿い単極。
flatnessによりclosed、従ってconstant-residue logarithmic。
detM成分ではβ=M−1ΘM−M−1dMの二重極を避けるためresθ1=resθ2。
ゆえにη=θ1−θ2の極はh成分だけで、Eulerη=d1−d2。
h constantのflags03/12ではη=0、d1≠d2に矛盾する。
nonlogなβをlogと仮定した証明ではない。
原証明と独立receiptは `PLANE_M2_ROOT_LINE_CLASS_OBSTRUCTION.md`、
`PLANE_M2_ROOT_LINE_CLASS_SEPARATE_REVIEW.md`。

## 3. flag11の全正常形は空又はMPS

hはline、detMはreduced quadratic。
irreducible conicならmarked root rowsをmixせずM=[[x,y],[z,−x]]へ戻せる。
Q1,Q2はrank四quadric threefoldでQ1∩Q2=S+T、Tはintegral cubic scroll。
一般ruling P1の留数和を使うと、hがどちらか一vertexを避ける全positionはinverse2formがdecomposableとなり矛盾。
両vertexを通る唯一位置h=xでは、全P0+P1+P2に対するweights(1,0,0,N,N)、N≥2がactualSNC curveを与える。
全higher sourceがMPSへ入る。

detMが二直線なら、全matrix imageの二/三次元を保つ。
image三はlegitimate row scales/normal右GL2で[[x,y],[x,z]]、image二はUdiag(x,y)VでleftU全parameterを保持する。
η=0より両場合とも一つのoff-diagonal connection entryは零。
rawPfのnormal u² coefficientは必ず3f k20だが、target係数はh x²のnonzero倍でfでは割れず矛盾。
全split-det sourceが空。

全proof/独立receiptは
`PLANE_M2_LINEAR_ROOTS_CONIC_MPS.md`、`PLANE_M2_LINEAR_ROOTS_CONIC_SEPARATE_REVIEW.md`、
`PLANE_M2_LINEAR_ROOTS_SPLIT_EXCLUSION.md`、`PLANE_M2_LINEAR_ROOTS_SPLIT_SEPARATE_REVIEW.md`。

## 4. constant-root三flagsは一様actualSNC又はnilpotent-curl空性

flags00/01/02は全primitive sourceを
M=[[1,0],[Aq,Bq]]、hBq=fへ戻す。
βはlower triangular、diagonal entriesはclosed logarithmic。
rawnormal bracketにはv²項が無い。

nonconcurrent triangleではbase weights零、normal weights(N+1,N)、N≥2。
transverse conic+line f=z(xy−z²)ではweights(0,0,1,N+2,N)、N≥4。
どちらも全P0非負・全P1/P2正で、actualDのweight零部分はκxyzuv、κ≠0。
したがって全sourceがactual GL5 curveでSNC pointへ達し、MPS Theorem1.4の全component scopeによりMPSへ入る。
period/resonance、全Aq、全higher coefficientsを固定したsliceではない。
必要な全weight boundsは `PLANE_M2_CONSTANT_ROOT_MPS_AND_NILPOTENT.md` §§1–3。

tangent conic+line又はconcurrent triangleでは、diagonal base fields X,Yはtracefreeの共通nilpotent spaceに入る。
Xf=Yf=0、X∧Y=0なのでactualG2のuv coefficient=3fκ、κ≠0。
fullP1のu-row bad coefficient v²Cに対し、projective Jacobiのnormalweight一係数とbasecurlから

\[
 (2\kappa I+4DY-DX)C=0
\]

となる。scalar+nilpotentでinvertibleなのでC=0。
fullu-rowがu倍となりactualDにu因子が生じ、既約性に反する。
off-diagonal Z、k20、他の全P1/P2はこのweight一係数へ寄与せず、全sourceで空性が従う。
全導出は `PLANE_M2_NILPOTENT_DIAGONAL_ROOTS_EXCLUSION.md`。

## 5. 全source表と未処理補集合

irreducible reduced cubic Rのrank二planeは既受理の全排除
`PLANE_IRREDUCIBLE_CUBIC_RANK_TWO_EXCLUSION.md` を使う。
残る四reducible Rと全六flagsは次の表で尽くされる。

| R | 00/01/02 | 03/12 | 11 |
|---|---|---|---|
| transverse conic+line | 全MPS | 空 | empty又は唯一marked位置の全MPS |
| tangent conic+line | 空 | 空 | 空 |
| nonconcurrent triangle | 全MPS | 空 | 空 |
| concurrent triangle | 空 | 空 | 空 |

これで指定全locusがMPS40の閉和へ入る。
SNC点の既知閉和所属だけをfamilyへ外挿せず、曲線を含む全ambient componentがSNC locusと交わるという定理を使う。
別higher scheme structure、全orbit分類、未知G全体の消滅を主張しない。

normal double quadratic Δ≡0、R非被約、m3、他support surfacesは明示残余。
分母はf及びdense normal-frame determinantを使う幾何に限られ、未知higher係数を反転しない。
全source正規化はactual linear inverses又はprimitive factorizationで戻せる。
rootcover/connection/source geometryとcurlは手証明、追加matrix/GB/jetの実行は零。
新constant-root/curlの正式独立receiptと復元済みconic/split receiptsを依存manifestへ固定して保存する。
