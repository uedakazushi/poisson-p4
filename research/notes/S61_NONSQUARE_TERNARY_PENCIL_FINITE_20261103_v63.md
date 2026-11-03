# Actual ternary pair の共通因子 locus は QQ 上で有限

2026-11-03。`A=X=W=0`（B自由）の実際の ternary pair Q0,Q3 を扱う。開は a*P*J!=0。今回証明したのは **二つの ternary quadratics が共通因子を持つ点の集合は零次元以下** というQQの有限性である。全域での空性やquad78は主張しない。

自然 core は `t=a/P, v=3e+20t, H=P+8t`。先行の元804項を使うsource replayから、実際のternary symmetric matrices M0,M3を再構成した。Q3の既知 row操作も保持する。共通lineがあればpencilの全conicがsingularになるので

`det(M0+s*M3)=c0+c1*s+c2*s²+c3*s³`

の全係数が零となる。逆は使わない。共通kernelを持つsingular pencilも含めた、この必要locus全体の有限性を以下で証明した。

既知unit t,P,Jと非零QQ定数だけでc_iを正規化した四多項式 C0,...,C3 の項数は28,56,91,131、全次数は4,6,8,9、t次数は2,3,4,4。C0を

`C0=L(H,v)*t²+M(H,v)*t+N(H,v)`

と書き、L=0の場合を別に保持した。

## 1. L!=0

各 C1,C2,C3 を C0 で pseudo-divide した。全ての literal identity

`L^k*Ci=Qi*C0+Ri`

を保存して検証した。Riのt次数は1以下。Riから非零QQ定数と、このchartで正当な t,P,J,L の因子だけを除いた。

最初の二つの t-elimination minors F(H,v),G(H,v) が互いに素であることを、巨大なresultantを展開せず次の二組の小さい恒等式で証明した。

- F,GのH係数から選んだQQ[v]多項式群に対する保存Bezout和が1。従って共通のvertical factorはない。
- v=0 で両式のH次数15,16が保たれ、specialized QQ[H]多項式に対する保存Bezout和が1。

共通factorにHの正次数があるなら、そのleading coefficientはF,Gのleading coefficientを割るため、両leading coefficientが非零のv=0でもH次数が保たれる。これはspecialized和1に矛盾する。よって全ての共通factorが除かれ、V(F,G)は有限。L!=0ではC0は非零の二次t多項式なので、その各(H,v)上のtも高々2点。したがってこのchartの必要locusは有限である。

## 2. L=0, M!=0

この場合、実際の C0 は `M*t+N` となる。このlinear polynomialに対するC1のpseudo-division identityを保存した。そのscalar remainder R(H,v) からMとJの因子だけをchart上で除いた。LとRについても上と同じ content和1 とv=0のspecialized和1を検証し、互いに素である。従って(H,v)は有限、M!=0なのでtは一意。このchartも有限である。

## 3. L=M=0

C0=0ならN=0も必要である。L,MのH-resultantおよびL,NのH-resultantを固定Sylvester determinantとして検証した。二つの一変数resultantに対する保存QQ[v]和は

`U*R_LM+V*R_LN=(v+3)(3v+4)/3`

となる。従って v=-3 またはv=-4/3のみ。各vでL,M,NのQQ[H] Bezoutを計算すると、それぞれ H-8、H+14/3 を得るので、H=8 またはH=-14/3となる。

この二点のそれぞれで、C1,C2,C3をtの多項式として評価した。三多項式に対する保存QQ[t] Bezout和が1となる。従ってこのleading-zero caseには候補がない。未知LやMを全域で逆にしていない。

## 範囲と再現

以上が全leading-coefficient caseを尽くす。したがってsingular pencilの必要locusはQQ上で零次元以下、その部分集合であるactual ternary共通因子locusも零次元以下である。proper finite Z_(p) model、全点の空性、必要locus上のactualquad78は未構成である。

この結果をquad78へ用いるには、訂正済みの別criterionにある **vertexbinaryrank2またはscrollbaseのgeometricintegrality** などの追加条件が必要である。旧omit-B²のconditional inferenceは別erratumでwithdrawnされており、本定理で復活させない。

探索計算は一つのbounded純QQ pseudo-remainder passと小さいcase completionで合計約8秒、最大expanded polynomialは16,402bytesで、指定1MBを下回った。GB・finitefield・追加primeを使わない。再現は既存 Cloud activation 後 `python code/s61_nonsquare_ternary_pencil_finite_20261103_v63.py`。元source、全pencil係数、pseudo identities、二組のcoprimality certificates、両leading-zero branchesを保存恒等式だけで検証する。
