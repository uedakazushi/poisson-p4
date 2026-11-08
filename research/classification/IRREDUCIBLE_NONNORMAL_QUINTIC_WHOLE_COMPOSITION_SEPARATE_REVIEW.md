# 全rank零planeと全既約非正規quintic：独立合成査読

2026-11-08、別作者 plane_collision_review。**純証明の合成査読 PASS**。
root の `PLANE_RANKZERO_CONDUCTOR_WHOLE_POINTWISE_COVERAGE.md` と
`IRREDUCIBLE_NONNORMAL_QUINTIC_WHOLE_POINTWISE_COVERAGE.md` を、全source分岐・適用仮定・
既受理の強いpointwise帰結に照らして独立に確認した。
新しい数学的算術、operator、kernel、minor、point、profile、GB、jetは実行していない。
旧G上界をpositive membershipへ読み替えない。新成分増分は零。

## 1. rank零plane入口の全first symbol

actual D は被約既約、Poissonのgeneric rankは四とする。
`REDUCED_DEGENERACY_ZERO_SURFACE_SEPARATE_REVIEW.md` の全m補題により、rank零planeの
conormal bracketは可換である。従って旧m二入口の全35/30/9 sourceを全normal位数へ戻せる。
旧m二だけで除外したleading零点を削除しないことが、この合成の前提である。

traceless co-Higgsの非零entriesはpointwise同一直線を張り、そのsaturationはline bundle O(k)⊂TP²。
非零係数sectionとH0(TP²(-2))=0によりk=0,1だけである。
O(1)のcoefficient map rank二以上は全co-Higgs定理へ、rank一はM=ell*Aなので
V=ell*Uというglobal vector fieldへ戻る。Vのdivisorial零を許すsupersourceを用いる。
O(0)も同じglobal Vへ入る。M全零/V全零はglobal V非零のunit chartに混ぜず、traceless零枝へ戻す。

global V非零のconstant sl2 matrixはnormalSS／normalnilの二型に尽くされる。
全weight零JacobiはSSでK=kzw,[T,V]=0、nilでK=az²+bzw,[T,V]=-bV。
SSおよびnil b零は可換なので、det(E,T,V)非零／恒等零を区別出来る。
前者はnormal位数二。後者ではconstant dependent T=tVと、独立だがpointwise dependentなcommon-imageを分ける。
全可換二field分類の他型Cartan/mixed/regular nil/common-kernelでは、独立ならこのdeterminantは非零。
constant dependenceとpointwise dependenceを同じsourceと扱わない。
normalnil b非零は可換分類へ入れず、正しい全四型・全rho・全high-m theoremへ戻す。

## 2. 全rank零plane表の強い帰結

| 全source | 合成に使うpointwise帰結 | G-onlyとの相違／全補集合 |
|---|---|---|
| O(1), rankM二以上 | rankM二は空、rankM三は既知閉和又は空 | 全kernel/nonkernel・Borel/ND・全tau・全m。正の全nil壁はactual第二rank二plane/m三へ |
| normalSS、独立T,V、Cartan | 全MPS閉和 | 全periodでSNC endpoint、全actual Cstar曲線とSNCを通る全component theorem |
| 同mixed三族／regularnil | 全空 | QQ2全壁とSS_MIX全三parameter、full badC/badD obstruction。独立別作者receiptを使用 |
| 同NILSQUARE/common-kernel等low21 | essential三八成分又は空 | 全generic/root、全七compatibility branch、FULL9 particular＋free三係数 |
| normalnil b非零 | 全空 | 全rho・八壁・唯一high-m rho=-1/3・残rho=-2/3の全actual shearとP2一意性 |
| normalnil b零、独立、det非零 | essential三八成分又は空 | mixed/Cartan全277階層、274必要核と全graph除外。残low21も全点処理 |
| 可換独立、det零/common-image | rank二以下 | 全四source、K零二源、K非零SS/nilの全shearと残N、全inhomogeneous P2 |
| dependent T=tV | essential三八成分又は空 | 全六source、H二源、全Ore八module、全central/zero補集合 |
| traceless零、T非零 | 全空 | 全15Jordan×K source、全20根、QQ2全allowed線／57根とoutside-domain復帰 |
| T零、K非零 | rank二以下 | K=zw/z²の両全型、平方型の残三方向も保持 |
| weight零symbol全零 | 全T36 | 全39実一軌道、全CE cocyclesの有限shear逆、今回のactual T36曲線と別作者review |

CartanのSNCはgeneric eigenvalueだけの主張ではない。
det(E,T,V)非零ならP0のactualDは非零定数×abc zwであり、periodの他の零を含めてlog-symplectic SNC。
全元Poisson componentがactual normal rescalingの極限P0を含むので、MPSの全SNC-meeting theoremを適用出来る。
「既知成分への一endpoint」だけから全元の所属を推測する議論ではない。

normalnil b零の混合/Cartan階層はselected-row必要kernelだけでbadC零を証明する。
完全kernelと新たに仮定せず、旧nonzero minor開上の正則な一意解と全真のwallを戻す。
人工的にdenominatorを払ったrank零点は元graph domain外／dependent lineであり、未知全kernelとして使わない。

各shear-exact証明ではfinite actual mapがP2九変数をtriangularに動かすことを保持する。
P1を消してからforcing零／D2単射を用い、homogeneous kernelだけをfull affine P2と取り違えない。
実解が存在する点だけでのP2一意性とactual stabilizerによる因子証明も、全parameterを含む正しい必要方向である。
元二乗係数からcomplex reduced supportでparameter零を得る場合も、全schemeの同型を主張しない。

したがって全rank零plane・全normal位数の各actual点が既知有限閉和又は矛盾へ入り、真の未被覆flagは無い。
旧G用8分岐の順序を用いるが、表ではtraceless零内の三補集合を明示的に分けている。

## 3. 全既約非正規Dに曲面支持があることとdegree上界

非正規integral hypersurface DはS2なので、R1の失敗を表す正規化商にintegral二次元支持Sが存在する。
一般二hyperplane切断を使う旧全入口から

```
sum(delta_i*deg S_i)=6-g,  delta_i>0, 0<=g<=5
```

を得る。従って各deg S_iは1..6である。
これは未知component一般元だけのdegree boundでも、g=0/1を全体へ仮定したものでもない。
一般平面曲線のnormalization genusと長さの式として、任意の元Dに適用する。

Hamiltonian derivationはchar0でnormalizationへ延長し、conductorを保つため各reduced integral支持はPoisson invariant。
D上でambient rank四は不可能、二次元支持上のgeneric rankは零又は二である。
少なくとも一つのSを選べば十分で、全支持を同時normal formに直すことは要求しない。

Sがdegree一ならplaneで、§1–2又は既rank二plane全位数合成へ入る。
Sがnonplaneなら旧全nonplane theoremを用いる。その全hypothesisを確認した：

* spanがhyperplaneでdegree>=3ならprimeCI I=(h,q_d)のsymbolic二乗はI²、degree五からh因子が生じ空。
  singular/normal性は要らず、degree二へこの比較を流用しない。
* degree二はsmooth quadric／quadric cone、両rankの全source。
  rank二coneのV1欠落を戻したFu-onlyを含む**V2とその別reviewだけ**を使用する。
* degree三はsmooth scroll／twisted-cubic cone、両rankの全四source。
  scroll rank零の撤回四sourceを使わず、正しい全L6＋Segre五sourceを使用する。
* degree四は全CI22／非CI quartic／projected Veronese／全centre。
  PVはG上界だけでなく全15 SL2/Borel点別空性、CI22はrank零の六源とrank二の十二源を区別する。
* degree五はlength式によりdelta=1,g=0/1。genus零で残るplane支持も保持し、
  rank零全113/122 projectionはfull Quot/部分導手変形からT40。
  ellipticの非ample／unstable境界も全sourceを戻す。rank二g0/ellipticは空、canonical-index二nodeは実T37、cusp等は空。
* degree六はdelta=1,g=0、全113/122/014/023/005 normalizationと全projection。
  両rankの全排除にsmall resolution/全split/非split/cuspを保持する。
  rank二023の旧proper-locus G上界ではなく、二つの実曲線による全residue矛盾を使う。

従ってnonplaneの全actual点はT37/T38/T39/T40又は空、all-planeも両rankの全plane定理で尽くされる。
normalization-flat、clean、ordinary crossing、reduced conductor又は単一支持という追加仮定は無い。

## 4. 依存の非循環性と査読の範囲

aggregate名だけでは相互参照に見える箇所を下位proofへ展開した。
coHiggs rankM三の正のsourceは、独立に完成した**rank二plane m三全点定理**へ送る。
nonplane A25は独立のconcurrent m二/T37実曲線を用いる。
plane tangent m二のt=-3は独立のdegree五rank零/T40 Quot theoremを用いる。
そのQuot theoremは全plane被覆又はnonplane全被覆を仮定しない。
これらは共通の閉MPS source／essential三定理へ向かう依存であり、aggregate同士を循環して証明していない。

一部入力（normalSS、scalar、zero-symbol等）のproducerは本査読者である。
それらの新数学をこの合成receiptで自己再認証せず、quartic_review、root、exceptional_coverageの
正式な別作者receiptを明示的入力とする。今回独立なのはrootによる全分岐の合成と仮定の適用である。

全既約非正規actual quinticの点別既知閉和への包含を受理する。新成分増分は零。
normal D、可約／非被約D、rank二以下のtensor、全非被約scheme構造、総成分上界、分類主定理の完成は射程外。
全evidence bindingとpure-proof gateはPASSだが、それ自体を数学の証明とは呼ばない。
元二notes／旧reports／主稿／global files／Gitは変更していない。
