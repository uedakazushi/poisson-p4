# B18 actual Dに一次で見えるnilpotentが存在する：DVR union補題

2026-11-03。新しい概念的証明。CAS・高次jet・既存producer再実行は用いない。
結論は指定P0のactual Dについて、nilradicalが極大idealの二乗に含まれないことである。
被約交叉の正確なembedding dimension、具体的なnil coordinate、nilpotence次数は決定しない。

## 1. 一般のfiber-product補題

Aを剰余体CのNoetherian局所環、B=C[[t]]をDVRとする。
共通の非零Artin局所商C0へ、A→C0とB→C0がともに全射であるとする。
Bのkernelは(t^k)、未知の有限k≥1であり、F=A times_C0 Bと置く。

F→AもF→Bも全射である。例えば任意のb∈BのC0像をAへliftすれば、そのpairはFに入る。
F→Aのkernelは

    K={(0,t^k b):b∈B}。

これはB-moduleとしてだけでなく、**F-idealとして**κ=(0,t^k)の一元生成である。
実際、任意のbに対してF→Bの全射から(a,b)∈Fを選べば、κ(a,b)=(0,t^k b)となる。
projectionの全射性を省略してprincipal idealと断定したわけではない。

Aの極大idealのe個の最小生成元を選び、そのF内のliftをx1,...,xeとする。
任意のm_F元のA像をこれらで表し、その係数をFへliftすると、残りはKに入る。
従って **m_F=(x1,...,xe,κ)** と厳密に成立し、

    edim F≤edim A+1。

κがm_F²に入る可能性はこの上界を弱くしない。generatorのliftをm_Fに取れるのは両射がlocalで、
K⊂m_Fだからである。完備生成元のliftだけから不正なtopological generationを推測したものでもない。
この補題は全k≥1に対するもので、有限個のkを計算した外挿ではない。

## 2. actual reduced sliceへの適用

受理済みのscheme-level orbit productから、

    D_hat ≅ S[[g19]],  edim S=20。

同じproductがN、Rとそれらのscheme unionにも同時に使える。
個別reduced Nの完成環はdomainでcone edim37なので、N slice A_Nはedim18。
Rはopen B18 orbitでregular20である。
このregular性は、既受理のgeneric B18におけるproper-source immersion/finite-Nakayama証明と、
P0のorbitがB18に稠密なopen orbitであることから得られる。
regular locusはG不変なので、そのgeneric証明はこのopen orbit全体に適用できる。
従ってR sliceはregular一次元のB=C[[t]]となる。

被約交叉の二つのminimal componentはNとR、全被約R∩NはB18。
G19を除いたR curveはN sliceに含まれず、両者は原点で会う。
よって、その**scheme intersection**の共通商C0はBの非零proper quotientで、(t^k)、有限k≥1となる。
このkは未計算のままでよい。
二つのcomponent idealのintersectionからの標準exact sequenceにより

    S_red ≅ A_N times_C0 C[[t]]。

単なる被約台の集合論的gluingではない。C0はN∩Rのscheme quotientであり、別のconductor商と
同一視していない。一般補題から

    edim S_red≤18+1=19。

## 3. 非零なlinear nil class

S→S_redのkernelをN_S=Nil(S)とすると、cotangentへの全射のkernelは正確に

    (N_S+m_S²)/m_S²。

actual edim S=20、edim S_red≤19なので、このkernelのC次元は少なくとも1。
従って **Nil(S) is not contained in m_S²**。
power-series十九因子を戻してもcotangentに加わるのはgroup方向だけなので、
**Nil(D_hat) is not contained in m_D_hat²** である。
同値にactual reduced completed Dのedimは高々38、actual Dは39であり、その差の少なくとも一方向は
一次で見えるnilpotentから来る。
edim S_red=19またはedim D_red=38という等号は結論しない。
T_RがT_Nに入るかどうかだけでこの結果を推論していない。

## 4. completionとreduced性の整合

元のlocal ringは有限型C-schemeから得られるexcellent局所環。
そのreduced quotientもexcellentであり、完備化はreducedである。
nilradical N0はNoetherian環でnilpotentなので、completionのexactnessから

    Nil(D_hat)=N0*D_hat,  (D_red)^hat=(D_hat)_red。

既知のreduced N/Rも同様に完備化できる。flatnessにより二component idealのintersectionは保存され、
上のfiber productは実際のreduced completed ringを与える。
cotangent spaceはcompletionで変わらないので、非零linear nil classは元のactual local ringにも降りる。
未代数化のformal nilpotentだけの存在を主張しているのではない。

この結果は、未決としていたnilradicalのlinear visibilityに対する新しい**存在**定理。
full原点D環、nilのexplicit representative、linear nil classの正確な数、D_redの正確なedim、
nilpotence・全primary構造・exact depthは依然別問題である。
independent conceptual reportは同名flat JSONへ保存する。
