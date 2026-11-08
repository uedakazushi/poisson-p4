# 複素射影四次元空間のPoisson構造：被約既約成分は正確に147

2026-11-08、integrator `/root`。Cloud開始commit4c3f623fb698f9adc6077015b25e2325b11072af、数学baseline a00a541b835b2bff18bafcba545ecbd4bae7175b。この再開区間の開始HEADは2b7c071f168d0709bb28e1cb436ec137d91b4200。元full150/350/S61全134と全actualPf、正逆、全補集合を保持した有限閉被覆と、各全成分性・相違証明を合成する。新しい五23成分を採録、最後の22候補は既知23成分の境界へ入り、新成分には数えない。

## 1. 最後のconstantTRIPLEの全actual Poisson族

補助座標(T,x,y,s)の積Poisson括弧{x,y}=1,{T,s}=sから出発する。
H=x³+3xy-3y/2,k=xy/3+y²/2-x4/12,A=e/12,B=e(15-e)/72,M=A x²y+B y²。

```
r=s+eT-H,
t=(1-ex/3)T-xs/3-k+M,
T=t+x(r+H)/3+k-M,
s=r+H-eT.
```

全eで両合成はidentity、全Jacobian determinant一。非線形な補助写像はPoisson族の構成にのみ用い、PGL帰属の代わりにはしない。
元のphysical全六括弧は

```
Pxy=1, Pxr=-3x+3/2, Pyr=3x²+3y,
Ptx=x/3+(1-2B)y-Ax²,
Pty=(-1+2A)xy+y/6-r/3,
Ptr=-et+(1-e/3)xr+(-3+3e/4+e²/8-e³/72)y²
    +(5e/12)xy-5y/4+r/2.
```

補助chainruleの高次x4/x²yは- e/4+3A=0、-e+eA-3A+6B=0で相殺し、全六がquadratic以下。e0は強制TRIPLEそのもの。signedPfはs=r+H-e[t+x(r+H)/3+k-M]、div=(-1+ex/3)∂t-e∂r。eunitではt係数-e、x4係数-e/4が非零で、fullactualDは単純h×既約quartic。multipleH/PfQ0のrankone17定理をこのgeneric simpleHへ無断適用しない。

## 2. genericfiberの既知全Ga9源へのliteral actual正逆

eunitで真のaffineGL4を

```
(X,Y,Z,w)=(t-r/e, ex/3, -r/e, 3y/e), det=1/e,
(t,x,y,r)=(X-Z,3Y/e, e w/3,-eZ)
```

と取る。modular場はY∂X+∂Z、P_Yw=1。u=X-YZでpotentialはΦ=u+φとなる。既受理GA9_CONSTANT_COMPONENTの完全十源parameterに

```
r_src=81/e4, b=3/e, alpha=-1/2,
p=q=D=G=k_src=0, A_src=-27/e4, H_src=1/2,
B_src=b-1/4, E_src=-b-1/6,
F2_src=e²(e-3)(e-12)/648
```

を代入すると、全六が一致する。とくにVy=1/2-bY,Vw=r_srcY²+bw、

```
P_YZ=bY-1/2, P_Yw=1, P_Zw=r_srcY²+bw,
P_XY=1/2+(1/3-b)Y-Y²/4+2F2_src*w,
P_Xw=Z+r_srcY²+(1/2-b)Yw+(b+1/6)w,
P_XZ=X-Z/2+(b-1)YZ+A2,
A2=F2_src(1+2b)w²-5Yw/12+5w/12.
```

r_srcは全eunitで非零。e3/e12のF2零等も元全十源で除外されない。これはpotentialの外形だけからの推測でなく、全源のactual座標正逆と全括弧へのliteral代入である。
従ってgeneric全fiberは既23成分R4_GA9_CONSTANT_HYPERPLANE_QUARTIC_23に入る。全physical係数はe多項式、全parameterでprojective係数非零（原点r1のPf1）なので、閉性によりTRIPLEも同じ成分へ。PGL不変性と受理済み全四source actualorbit正逆を戻せば、元highesttriple枝の全sourceが既知閉和へ入る。

独立root再生 `CLASSIFICATION_CONSTANT_MODULAR_TRIPLE_GA9_FULL_FAMILY_ROOT_REVIEW_20261108.json` は元積Poissonから全chainrule、両逆/det、全150canonical係数、全350原JacobiとstrictJacobi、全126actualD、e0、genericGa9frame/Phi_uを確認してPASS（0.315305秒、RSS50072KiB）。別作者は全sourceの適用とliteralc1正規化を純査読した。算術再現とpure読解を区別する。

## 3. 全被約支持の有限閉被覆を閉じる

旧global149全case監査は、全rank四点を既142及び七既約閉像へ送る：Ga9linearw、plane二H+cubic三源、constantDOUBLE/TRIPLE、unimodP*。全actual五次素因子分割、全Jordan/coker、全lowerfiber、全zero/unitと非飽和thetaを保持している。
P*→真T19で全5Hを既知和へ、七→六。Ga9/plane三/DOUBLEの五23閉像は、wholeT36/T38 properflags、wholeT37単純因子degree≤5条件、全prior142の次元/closedfactor/接触/actualgraph監査及び候補相互非包含から全成分と証明した。既142と相異なる五つなのでknown147となる。
最後に§1–2が残TRIPLEを既Ga9constantへ送り、全rank四点がknown147へ入る。known147は有限閉和だからrank四全閉包も入る。rank二全支持 theoremの「六rank二成分又はrank四閉包」のclosure項へ**この後で**代入して、全rank二点も同じ和へ。rank四被覆より先にrank二closure項を循環使用しない。

従って複素P4のPoisson parameter space（PGL商前）の被約支持について

```
P_red = union of the 147 individually catalogued closed irreducible components.
```

全147の個別成分性と相互相違は既監査＋今回五component証明による。任意の他全成分は有限和の一項へ入り極大性でその項に等しいので、追加成分は無い。rank四141、rank二6、全成分数は**正確に147**。

次元別の数はrank四で22:18、23:89、24:12、25:20、28:2。rank二で22:2、28:2、36:1、39:1。計141+6=147。一般元の全actual construction・退化因子幾何・次元は各component recordと同一主原稿へ記載する。

## 4. 完了の射程と再現

主分類の全支持・成分性・相違・一般元の義務は完了。未処理primarysource、unknowncandidate、分母の未確認補集合は無い。各元定理の仮定とexcludedchartはそのまま保持し、simpleHとmultipleHを取り違えない。
全orbit同型類、全schemeの非被約厚み、正確nil指数/H0/T2等のsecondary問題を解いたとは主張しない。新五成分のsmoothness/局所唯一性も追加仮定にしない。plancheck/gate/testのPASSを数学証明そのものと混同しない。
最小依存の再現入口は本note、全case旧globalcover、五component/root査読、P*/TRIPLE同族全係数再生、各独立source/certificate manifest。旧producer/実行済みcode/report/job/manifestは不変。newliteralcallerはinline全sourceを保存し、既reportsを上書きせず別work出力先だけを変更して再生できる。環境は既CloudのM2 1.24.11/SymPy/Fraction/LuaLaTeXを再利用する。
本receiptは数学の完了を記録する。主稿最終組版・可視独立確認・通常main commit/pushの実状態はそれぞれ専用receiptとGitで別に確認する。
