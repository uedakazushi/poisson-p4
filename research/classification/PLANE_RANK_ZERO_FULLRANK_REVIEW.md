# rank零plane・二重因子・O(1)係数rank三の全像上界：独立レビュー

2026年11月7日。主義務G。**全m_D=2、k=1、rankM=3枝の射影Poisson像はdimension≤23**
という主張を、全源の幾何的還元と別作者Fraction算術の両面で受理する。
rankM≤2、k=0、multiplicity≥3、rank二restrictionの枝はこの結果に含めない。
全E帰属又は全Poisson成分分類の完了も主張しない。

独立checkerは `code/classification_plane_rank_zero_fullrank_review.py`、
証明書は `reports/CLASSIFICATION_PLANE_RANK_ZERO_FULLRANK_REVIEW_20261107.json`、
有限契約は `jobs/plane_rank_zero_fullrank_review.json`。
producer、producer helper、SymPy、CAS、NumPy、legacy算術をimportしない。
既存の別作者標準Fraction疎多項式演算と一変数QQ(τ)演算だけを再利用した。
gate後、2.457秒、最大RSS22,164KiBで全検査PASS。
最初の実行はPythonの条件式の空白不足によるsyntax errorで算術開始前に停止した。
その一箇所の表記を修正して契約を再gateし、全算術を実行した。

## 1. planeへのFrobenius入口と全係数を保持するEuler商

`RANK_ZERO_SURFACE_GENERIC_FIRST_JET.md` §§2–3を独立に読解した。
同noteの「Sがplaneでない」仮定は§1でm_D=2を導くためだけに使われる。
今回はm_D=2を明示仮定するため、§§2–3の代数をplaneへそのまま適用できる。

一般点の四次元cotangent Lie algebra Lは二次元conormal ideal Nを持ち、商は可換である。
πがplaneで零、Dのnormal二次項が非零なので、linearized PoissonのPfaffianは非零、LはFrobeniusである。
Nがnonabelianならaff1であり、全derivationがinner、centerが零である。
商の二liftからinner derivationを除くと両者はNと可換になり、その相互bracketもcenter零から零。
従ってL=aff1⊕C²、rank≤2となって矛盾する。よってconormalはabelianである。
一般点で零のbundle bracketは全plane上で零となる。
first-jet noteの「非Frobenius補集合」はこのm_D=2内の未処理点を意味しない。

planeを(z=w=0)とすると、全projective零性はraw90係数、Euler冗長十でdimension80。
conormal可換性でnormal-normalのnormal次数一の六係数を除き、weight零35、一30、二9となる。
weight零はH⁰(TP²)⊗End(C²)の32とvertical二次Kの三係数である。
weight一ではbase-normalのz²∂z、w²∂wの六係数を、Euler項で一意に除く。
これがcheckerで再構成した全30 sourceであり、weight二の全九係数も保持する。

## 2. 全co-Higgs線とLax還元の必要十分性

Jacobiのbase二・normal一部分はtraceless三entryのpairwise wedgeを零にする。
その共通imageのsaturationはTP²のline bundle O(k)である。
非零global sectionからk≥0、Euler列によるH⁰(TP²(−2))=0からk≤1。
今回扱うk=1ではinclusionはH⁰(TP²(−1))=C³の非零vector U、すなわちpoint p∈P²で指定される。
選んだpのchartでU=∂x0とし、全sourceを

    Π0=T∧E_n+U∧M(x)n+K∂z∧∂w

と書く。Tはtrace-free matrixの全八係数、Mはlinear sl2の全九係数、K=q_Sは全三係数を持つ。
checkerはこれら二十係数を記号のまま保持し、strict Jacobiの全vector係数が

    q_M[T,U]+U(V_MK+q_T(M)+det(Mn,M0n))

に厳密一致することをFraction疎多項式で再構成した。
他の七strict Jacobi成分も恒等零である。
これは生産者のLax式を入力として演算した検査ではない。

[T,U]がUと独立なら上式は一定の二次元vector空間に入り、
全xでEuler_baseに比例するには恒等零でなければならない。
独立二係数を比較するとM=0となってrankM=3に反する。
従ってTはpを保存し、[T,U]=−t00U。
二×二traceless恒等式を代入して

    T(M)−t00M+[R,M]=0, R=S+M0/2

を得る。逆代入も同じJacobi恒等式を零にするので必要十分である。

rankM=3ではMの逆が存在し、T=t00I−M⁻¹ad_RM。
trace零からt00=0、x0係数から[R,M0]=0。
M(p)≠0はMの単射性とp≠0から従い、そのsl2 centralizerは一次元なのでR=τM0。
baseの逆行列M⁻¹、normal GL2の逆及び非零tensor scalarにより、
非零M0のsemisimple・nilpotent二軌道へ戻せる。
semisimpleの非零固有値の規格化は有限代数cover上で行えば十分で、次元を増減しない。
higher-weight全係数もこれらblock変換で可逆に移り、隠れたmodulusや捨てたfiberはない。

標準M=[[b,c],[a,−b]]に対する二模型のT,Kは、ad_Hとad_Fを直接計算して
producerの表示と一致する。
旧記号c=0は現在のτ=0と同一であり、base座標c又はM=0を意味しない。

## 3. 全350行、全τ、真正shearの別算術

二模型のΠ0、独立30個のΠ1単項式、独立九個のΠ2単項式から、
Schouten括弧とEuler wedgeを別に実装して全350行の二行列を再構成した。
保存行列350×30及び350×9と全係数が一致する。
actual coordinate shear `x→x+zY+wZ` の全六微分も別に作り、
上記の真正Euler正規化後の30×6行列と一致する。
全D1行へのshear代入も零である。

QQ(τ)上の完全kernelを全行へ代入して零を確認し、kernel独立性も別行消去で確認した。
選択maximal minorの行列式はFraction有理関数による別Gaussian演算で計算した。
下界minorと独立kernelによる上界が一致するので、generic rankは推測ではない。

全minor積が、保存されたlinear factorsの積に非零定数倍で一致することを検査した。
従って未確認の非有理又は代数的例外根は残っていない。
ss十根、nil六根の全てについてQQ上のminor、完全kernel、shearを再検証した。
τ=0も算術上は保持し、m_D≥3への明示補集合として扱った。

両型のgeneric nullityはD1が六、D2が零。
τ≠0の危険な壁はssの1/2で(6,1)、nilの1で(6,2)、3/2で(8,0)、1/2で(7,1)であり、
他の壁も全て証明書の表に保持した。
この結果は全matrix上の有理計算であり、幾つかのτ標本を一般化したものではない。

## 4. 実D、τ零の除外、全parameterの復元

Euler収縮した五つの符号付きPfaffianを独立に作り、

    D_ss=−4acτ(az²−2bzw−cw²),
    D_nil=c²τ(az²−2bzw−cw²)

を全係数で確認した。
これはnormal次数二のleading Dそのものであり、higher-weight項は次数三以上にしか寄与しない。
τ=0ではleading Dが零となるためm_D=2ではない。
τ=0をemptyと呼んだり、位数三以上の分類を済ませたとは扱わない。

fixed planeのfirst symbolはp二、M九、τ一でssが十二次元。
nil条件det M(p)=0はM可逆・p非零のopen上の一つの真のdiscriminant方程式なので十一。
τは独立なparameterであり、finite wallへ固定すれば各々一だけ下がる。
fixed Π1のΠ2条件はD2のinhomogeneous linear systemなので、
全higher Poisson fiberの次元はdimkerD1+dimkerD2以下。
未解決の高weight非線形式はこのsupersourceを狭めるだけである。

plane選択六、tensor projectivization一を戻すと、generic上界はss23、nil22。
fixed wallではssが16+d、nilが15+dで、全τ≠0の壁でも二十三以下となる。
全firstsymbolの可逆normalizationでfiber dimensionを移送しているため、
一点のslice次元を全像へ外挿していない。

genericではD1 kernelがactual shear六と一致し、D2は単射である。
shearはnormal weightを一つ上げる実unipotent linear transformationで、
二次bivectorの最大weightが二だから有限conjugationによりΠ1を消去できる。
変換後のweight二JacobiからΠ2=0。
上のD因子によりこのgeneric枝は既約五次を持たない。
壁上ではこの強い空性を使わず、全parameterの弱いdimension上界だけを使う。

以上により、明示されたm2/k1/rank3枝の**全像≤23**を受理する。
成分数、個別成分所属、rankM≤2又は他の補集合について結論を追加しない。
