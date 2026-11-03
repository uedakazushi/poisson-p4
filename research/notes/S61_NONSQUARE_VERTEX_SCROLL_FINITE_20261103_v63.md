# Scroll の integral 開と vertex 比例 locus の有限例外

2026-11-03。actual raw four quadrics の source、座標 `(A,X,W,B,S,U)`、`c=6075` は凍結済み binary reduction と whole-pair packet に従う。基礎開は `a*P*J!=0`、自然 core は `t=a/P, v=3e+20t, H=P+8t` とする。ここで t と P は unit である。

このノートは旧 omit-B² 条件からの誤った syzygy 推論を使わない。旧条件は明示的 erratum により撤回されており、以下の帰結は別の **corrected conditional** の `G1∧G2∧(G3 または integral C)` を用いる。

## 1. 実際の scroll と rank-drop resultant

厳密な表示は

`Q1=AB/c+W*C1`, `Q2prime=40X*B/c+W*C2`, `Dquad=40X*C1-A*C2`。

`C1=alpha*A+beta*X+gamma*W+delta*B`,
`C2=epsilon*A+zeta*X+eta*W+theta*B`

で、actual coefficient は `beta=32*a³*J²/243` という既知 unit を持つ。三式は `[[A/c,40X/c,W],[C1,C2,-B]]` の maximal minors（符号・定数倍）である。既に Q1,Q2prime は互いに素と証明され、minor ideal は height2。Hilbert–Burch の最小分解から、base `C⊂P³_(A,X,W,B)` は CM、純一次元、次数3であり、その座標環に一次式は含まれない。

第一行が非零の点では第二行との比を lambda として `B=-lambda*W`。残りの二条件は

`M(lambda)*(A,X,W)^T=0`,

`M(lambda)=[[alpha-lambda/c,beta,gamma-delta*lambda],`
`           [epsilon,zeta-40lambda/c,eta-theta*lambda]]`。

beta が unit なので、この 2×3 行列の rank<2 は

`r1=beta*epsilon-(zeta-40lambda/c)*(alpha-lambda/c)=0`,
`r2=beta*(eta-theta*lambda)-(zeta-40lambda/c)*(gamma-delta*lambda)=0`

と同値である。literal actual coefficients から、分母を既知 P²だけで払った二つの quadratic の固定4×4 Sylvester determinant を照合した。結果は

`Res_lambda(P²*r1,P²*r2)`
` = -67108864*a^9*J^6*P²*Eold /4749088069944112812890625`。

自然 core では `Eold/P³=E5` で、

`E5=2*(v+2)*(7v-11)^2-t*(55H-66v+18)^2`。

したがって E5!=0 上では全有限 lambda で rank M=2。E5 は未知 case equation であり、全開の unit として扱わない。

## 2. E5!=0 上の C は幾何学的 integral

行比を `[s:tau]∈P¹` として斉次化し、

`M(s,tau)=[[s*alpha-tau/c,s*beta,s*gamma-tau*delta],`
`          [s*epsilon,s*zeta-40tau/c,s*eta-tau*theta]]`

とする。s!=0 では上の resultant により rank2。s=0 では最初の二列の determinant が `40*tau²/c²` で非零だから、ここでも rank2である。

二行の cross product を `(Ac,Xc,Wc)` とする。これは共通零点を持たない三つの次数2形式であり、s=0 で `Wc=40*tau²/c²`。よって

`phi:P¹→P³, [s:tau]↦[s*Ac:s*Xc:s*Wc:-tau*Wc]`

は basepoint を持たない三次形式の morphism である。有限比では上の kernel 方程式から全 minor が零になる。無限比の像は `P0=[0:0:0:1]` で、そこでも第一行が零なので C に入る。

逆に C の点で第一行が非零なら比 lambda が一意で、その (A,X,W) は rank2行列 M の唯一の kernel line になる。従って phi の像に属する。第一行が零となる点は P0 のみであり、これも無限比の像である。したがって phi は C の underlying support を尽くし、一般点で一対一である。

basepoint がないので `phi^*O(1)=O_P¹(3)`。一般次数1から像は irreducible degree3 curve。Hilbert–Burch scheme C の次数も3だから、C は一般点で長さ1、従って generically reduced である。CM は S1、従って generically reduced と合わせて reduced。これで C は integral となる。この議論は任意の代数閉 characteristic-zero residue field で成り立つので、geometric integrality である。smoothness や smooth⇒integral の省略は用いない。

## 3. G3 failure と E5=0 の交わりは有限

G3 は vertex `A=X=W=B=0` 上の Q0,Q3 の binary coefficient matrix の rank2条件である。実際の pair monomial order は `[B²,B*S,S²,B*U,S*U,U²]` なので、使う列は正確に `[2,4,5]` である。B*S/B*U 列を用いない。

whole-pair source replay の literal minor normalizations から、G3 failure には二つの式

`alpha_i(H,v)*t+beta_i(H,v)=0 (i=0,2)`

が必要である。元 minor から除いたのは非零QQ定数と t の既知 unit のみ。特に f=(v+3)(3v+4) を割っていない。

`D0=55H-66v+18`, `N0=2(v+2)(7v-11)^2` とする。E5=0 は `D0²*t=N0`。

D0!=0 では t が一意となり、必要な二つの projection polynomials は

`Fi(H,v)=D0²*beta_i+N0*alpha_i`。

二式の項数は26,21、H次数3,2である。保存証明書は次を literal QQ identity として照合する。

1. 二式の H係数に対する QQ[v] Bezout 和は1。従って共通 vertical factor はない。
2. v=0 で H次数3,2が双方保たれ、specialized QQ[H] Bezout 和は1。

共通因子の H次数が正なら、その leading coefficient は両式の leading coefficient を割る。v=0 で両 leading coefficients が非零なので、その因子の H次数は保たれ、2に反する。従って二式は互いに素であり、その共同零点は有限。各点で t は一意だから、この chart の locus は有限である。

D0=0 では E5=0 から N0=0。従って `(v,H)=(-2,-30/11)` または `(11/7,120/77)` の二場合のみ。各点で上の二つの actual linear-t equations に対する保存 QQ[t] Bezout 和が1なので、双方とも空である。f=0 も第1 chart に含まれ、除外した未知因子はない。

したがって `V(E5)∩{rank vertex pair<2}` は aPJ 開で零次元以下である。これは QQ の有限性であり、全域の空性や finite-Z_(p) model は主張しない。

## 4. Raw quadratic-cofactor rank78 の有限欠損

corrected conditional packet における G1 は Q0,Q3 mod(B,W) の互いに素、G2 は Q0,Q3 mod(A,X,W) の互いに素である。

G1 failure は既証明の finite cover Sigma に含まれる。high-rank proportionality は aPJ 上で literal QQ identities により除外済み。rank2 branch は DB=0 と二多項式 N,M の共同零点に入り、monic e4 equation と monic P28 membership によりこの closed cover 自体が QQ 上で有限である。この有限性だけには F,B0,C の可逆性も、raw24の finite-Nakayama conclusion も不要である。

G2 failure は別の frozen ternary-pencil finite theorem により aPJ 上で有限。G1,G2が共に成り立ち、さらに今回の有限集合 `V(E5)∩G3-failure` を除くと、G3が成立するか、E5!=0 により C が integral である。従って **corrected conditional** が適用され、四 raw quadrics の quadratic-cofactor map `84×126` は rank78となる（六つのKoszul syzygiesだけ）。

以上から raw rank78 defect は **全 aPJ 開で零次元以下** である。有限集合の各点で rank78かどうかは未検証であり、全点下界とはしない。旧 omit-B² の撤回済み推論は使わない。

## 5. 元の六二次式の次数4への帰結

元の six-quadratic ideal を I、四式の ideal を J0、ell と Z,W の二式を ell*Z,ell*W とする。FULLU 上で actual ell chart が可逆である。raw quadratic-cofactor syzygies が六Koszulだけなら、`sum Qi*Ai2=ell*f3` を mod ell に取ると cofactor は constant Koszul のみ。それを元の四式間の Koszul で引いた後、各 Ai2 が ell で割れ、次数4での colon は必要な範囲で J0_3 に一致する。

従って I4 の mod ell image は rank78、kernel は

`ell*(J0_3+(Z,W)*S2)`

となる。既証明の original Z=W=0 five-variable product rank19 により quotient の次数3は `35-19=16`。七 light variables の次数3は84だから kernel rank は68。結局 original light quadratic-products の次数4 rank は `78+68=146` である。

したがって FULLU では有限集合を除き light rank146、H倍の独立六方向を加えた weighted quadratic-product rank152。全 core での upper146/152 は既存 packet に保存済みである。全点の lower146、full cubic57、残る actual quartic quotient 全射、横断長65、ambient除去は、この有限欠損定理からは結論しない。

## 再現と保存

`python code/s61_nonsquare_vertex_scroll_finite_20261103_v63.py <new-output>` は既存 Cloud tools を再利用し、actual804 source、four-coordinate coefficients、Sylvester determinant、E5変換、二つの正しい vertex minors、content/specialization Bezout と二つの D0=0 branches を照合する。GB、新しい有限体、大きな行列展開は行わない。幾何学的 coverage・CM/reduced の議論と corrected conditional の適用は本ノートの証明である。
