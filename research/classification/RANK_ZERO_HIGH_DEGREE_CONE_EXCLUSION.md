# degree五以上の全cone支持：rank零・既約quinticを排除

2026-11-07、主義務G。計算を使わない全source補題。
`PROJECTED_QUARTIC_CONE_RANK_ZERO_EXCLUSION.md` §2のEuler-weight導出を一般化する。

**命題。** S⊂P4をpoint vertex vを持つintegral・非退化cone surfaceとする。
S上で零となる任意のglobal projective bivector πについて、
base curve C⊂P3のquadratic idealが零ならrankπ≤2、一dimensionalならactualDは
そのunique quadricで割れる。
特にdegree S≥5なら、rank四Poissonかつ実Dが既約被約quinticとなることは無い。
πのJacobiを必要とせず、normality、base curveのsmoothness、linearly normal性、
arithmetically Cohen–Macaulay性又はconormal bracketを仮定しない。

## 全Euler sourceの証明

v=(0:0:0:0:1)、base hyperplaneの座標をx=(x0,x1,x2,x3)、vertex座標をcとする。
integral非退化curve C⊂P3について、Sは(x,c)、x∈affine cone(C)でparametrizeされる。
全global bivectorは、Euler商の前に二次係数のraw liftを持つ。
このliftを全て保持して

`Π=P(x)+cL(x)+c²A+∂c∧(V(x)+cW(x)+c²v0)`

と書く。P,Vは二次、L,Wは一次、A,v0は定数である。
Euler E=E_x+c∂cとのwedgeは

`E_x∧(P+cL+c²A)+∂c∧[cP+c²L+c³A−E_x∧V−cE_x∧W−c²E_x∧v0]`。

S上でprojective πが零なので、各c係数はC上で零である。
c³係数はA=0を与える。
c²係数は `L=E_x∧v0` をC上で与える。
これは一次係数の恒等式であり、CがP3をspanするのでambient全体に成立する。
従ってEuler gauge `E∧(c*v0)` を引けばLとv0を同時に除ける。

c係数は `P−E_x∧W` の全二次係数がI_C(2)に入ることを与える。
Euler gauge `E∧W` を引くことで、mixed cW項も同時に除ける。
従って全sourceは

`Π=Σ_j Q_j(x)*A_j+∂c∧V(x)`                         (1)

の代表を持つ。Q_jはI_C(2)のbasis、A_jは任意constant skew bivectorである。
残るc零係数は `E_x∧V|C=0`、すなわちVがC上でradialという正確な条件である。
Vをquadric-times-constantへ制限していないので、nonCM extra sectionも保持する。

I_C(2)=0なら(1)はdecomposableな `∂c∧V` だけとなりrank≤2。
I_C(2)=C*Qなら

`Π∧Π=Q*(Q*A∧A+2*A∧∂c∧V)`。

actual anticanonical DはE∧Π∧Πの五vector係数であり、従ってQで割れる。
これはEuler gaugeに依存しない実Dの因子である。
rank四ならDは非零なので、degree二Q因子は既約quinticと矛盾する。

## degree五以上での全網羅性

degC=degS≥5とする。
I_C(2)に独立な二quadricがあるなら、integral非退化Cを含むquadricはreducibleではない。
実際、linear factorの積がprime I_Cへ入ればその一因子がI_Cへ入り、Cがplaneに入る。
従って独立な二quadricはcommon factorを持たず、degree四のCI curveを切り出す。
それがdegree≥5のCを含むことはdegree比較に反する。
よってdimI_C(2)≤1であり、上の二caseで全て排除される。

この結論はdegree五・六の非normal cone及びnormal coneの双方を含む。
既済normal-degree-five認識又はnormal-sextic ruled還元の仮定を要しない。
non-cone支持とπ|Sのrank二枝は補集合として残る。
