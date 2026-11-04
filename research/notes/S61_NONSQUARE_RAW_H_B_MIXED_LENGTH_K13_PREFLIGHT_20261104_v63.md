# 共通monic基底のFitting長さ上界13とK13入力preflight

2026-11-04。混合normは別の独立review packetの受理範囲を参照する。本稿自身はそのreviewではなく、受理されたsource/PID補題を保存し、次の有限係数代数計算の入力だけを固定する。新しいnorm、e-Euclid、pivot逆、field Gaussian、点、素数、GBは実行しない。

## 1. 長さ上界は必要a値の単純な数え上げではない

R=QQ[a,a^-1]、C=R[e]/(h1)とする。二つのsource inputは全く同じe-monic16のh1と基底1,e,...,e15を持ち、Mk1とMm1はそれぞれk1とm1の乗法行列である。全成分に負のa冪はなく、整数clear Dold,Dnewは非零QQ定数だけ。従って

T=coker([Mk1|Mm1])=R[e]/(h1,k1,m1).

Fitt0(T)には両方の16次square subdeterminantが入る。保存された全QQ Bézoutとcontent/a-power除去により、primitive13次g13もFitt0(T)に入る。一方Tは既存長117の代数の商なのでQQ有限。PID RでTのorder polynomial d_TがFitt0を生成し、d_Tはg13を割る。従ってQQ-length(T)≤13。これは13個のa値に16を掛ける粗い評価ではない。

既存HはR[e]上d-monic4だからT[d]/(H)はT-free4で、明示coverのQQ長は高々52。ここでTをactual Bのcokernelと同一視しない。

## 2. actual bad schemeとの関係に必要なgate

actual Bbad上でこの三式が必要なのはD(G)に限る。G0は既存Q/L例外に保持する。

schemeとしても52上界を使うには、mixed係数のsource ideal所属と全unitを束縛する。有限recipeは次である。モニックHの四次代数でA(lambda)=s*(Mg+lambda*Mk)、W(lambda)=adj(A(lambda))*e0を取る。四本のsource monic除算をliftすると

det A(lambda)=s*(g+lambda*k)*W(e,lambda)-s*H*Q(e,lambda).

lambda1係数を取ればraw mixed coefficientは(H,g,k)に属する。保存unit-strip因子はFULLUの単元であり、さらにD(G)でだけGを除く。これらを用いて初めて明示coverからactual restricted bad schemeへの商写像を述べられる。全G0を含むactual Bbad schemeの長さ52は主張しない。今回mixed adjugateのscalar係数自体を展開したわけではない。

## 3. K13は既存factorを再利用し、単元を実証して計算する

g13は以前のOLD_A117_FACTORに保存されたdegree13 exponent1 factorと非零QQ倍で一致する。全係数を比較し、そのQQ単元を保存する。新しいfactorizationは行わない。旧native factorのirreducibility provenanceを引用するが、今回独立irreducibilityレビューを行ったとは呼ばない。

将来の計算はQQ[a]/g13の13成分vectorで行えば、irreducibilityを黙って仮定しなくても安全に実装できる。pivot p(a)を反転する前にxgcd(p,g13)=1とS*p+T*g13=1の全QQ係数を確認する。非単元なら両成分を保持して分割し、零pivotは飛ばす。a自体はg13 constant termが非零なので、a^-1=-(g1+g2*a+...+g13*a^12)/g0という明示identityを使える。

## 4. 小さいe入力と次の完了条件

既存の二つの乗法行列の第一列は、k1 mod h1とm1 mod h1そのものである。これを再利用すればe入力はh1の16次、二つの15次以下remainderで足り、25次k1からやり直す必要はない。今回aだけをg13でmonic reductionして、全係数を13成分以下へ正確に落とす。この操作の逆元はg13最高係数というQQ単元だけである。

次に予算判断を求めるのはgcd_e(h1,m1remainder)、続いてold k1remainderとのgcdという一回の有限代数計算。degree0なら単位Bézoutで明示Tが零となる。degree1ならeがK13上有理になるが、その先のmonic H4と全22 source条件の有限NF検査は別予算。高い次数なら有限quotientを保持する。

input数・実coefficient height・最大663 field coefficient slotsを保存する。Euclid中のheightが小さいことを未証明のまま保証せず、将来jobには60秒/2GB/32MiBのcheckpoint guardと単元証明を課す。今回そのjobは実行していない。

## 5. 優先する、field inverseを使わないmixed-minor判定

同じh1乗法なのでNkとNm1は可換。旧V=adj(Nk)e0は、独立受理済みNV=det(Nk)e0と非零detから厳密に定まる。Vから**全16成分共通のQQ contentとa冪だけ**を除き、Vnormalとする。各成分を別々にscaleしてはいけない。

W=Nm1*Vnormal mod g13を考える。可換性からNm1*adj(Nk)=adj(Nk)*Nm1。従ってcommon unitを戻したW_iは、Nkの列iをNm1の列0で置き換えたdeterminant、つまり[ Mk1 | Mm1 ]のmixed16-minorである。既知単元を除いてもFitt0(T)への所属は変わらない。

よってgcd(g13,W0,...,W15)=1ならFitt0(T)は単位でT=0となり、旧D(G)部分を除去できる。このcriterionはe-EuclidもK13のfield inverseも使わない。W全成分が零でもlength13やe有理性をまだ結論しない。

今回入力preflightでは、Vのcommon normalization、全16個のmod g13剰余、Nm1の全256成分のmod g13剰余だけを構成する。**Wの積は計算しない。** input factor13のirreducibilityを使わずとも、gcdと全QQ Bézoutが単元criterionを証明できる。

次の一回を提案するなら、16×16のunivariate多項式行列と16vectorの積・mod g13・degree13 gcd証明だけでよい。上界は43264個のscalar convolution積と2496個のreduction係数更新、Wの出力208係数。common numeric denominatorとprimitive整数g13によるpseudoremainderを使えば、未知parameterの逆元は不要である。実input heightとそれから導くWのbit/保存上界を別JSONに記録する。これは新計算の許可ではなく予算前flightであり、e-Euclidより先に評価する候補である。

## 今回の入力結果と次の一回の規模

preflightは1.9908秒で完了した。h1、二つのe-remainderは81/159/104項、e次数16/15/15、a次数8/12/9。最大QQ numerator/denominatorは241/230bitで、e方向の計算は行っていない。

旧Vの全成分共通contentは4433bitの整数、共通a冪は24。これらだけを除きmod g13すると、**全16成分が非零**で計208項になった。Nm1 mod g13は2529項。numeric共通分母はそれぞれ311/8359bit、整数化係数は634/9935bitである。独立なcolumn scalingは行っていない。

整数g13の最高係数Lを使った12段以下のpseudo-remainderから、Wの整数numerator上界12629bit、numeric denominator上界10638bit、RAW上界1488032 bytesを得る。degree13と非零W成分のxgcd certificateも、25次以下の固定Sylvester/cofactor heightから数MiBの範囲で提案可能である。ただし今回Wやxgcdは実行していない。全input artifactは約2.38MBである。

次に提案する一回は、この保存入力だけを使うNm1*Vnormalの全16成分、整数pseudo-remainderの全商・積恒等式、g13とのgcd/QQ Bézoutの確認である。60秒/2GB/32MiB、零Wまたは非単元gcdなら全成分を保持し、別点・素数・m2・field Euclidへ自動切替しない。
