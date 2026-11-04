# A0 normal block と疎 cofactor map の読出し

2026-11-04。既存 `S61_AFFINE_RESIDUAL42_20261103_v63` の同じ QQ 中心・線形129×42行列だけを読む。D1再計算は行わない。`S61_AFFINE_NORMAL_LEFTINVERSE_ASSESSMENT_20261103_v63` は既にD0を一般betaで排除し、同じmod32003中心では完全Macaulay計算からD1を排除している。後者をQQや一般betaの障害には昇格していない。その既存結果を重複した新しいnegative theoremとして数えない。

新しい小さいデータは、線形factor matrixの正確なsparsity、first-normal規則によるtriangular support、cofactor-mapの実際に非零なconstraint rowsと初期nonzero数、定数42 blockの正規化recipeである。全行列を再生成せず、coefficient solver・GB・新しいfield/pointは使わない。

変数は(z0,z1,z2,n0,…,n41)。Aのnormal column jでは、normal変数n_kの係数はk≥jにしか現れない。これは単項式の最初のnormalで因数化した保存規則そのものである。従ってD1 coefficient mapの名目shape45402×5934には構造的なzero rowsがある。全45変数のcofactor multipliers1,q0,…,q44を使う incidence を数えるが、実際の大行列は組み立てない。dense fmpqのメモリ見積りはentry本体だけであり、rational numerator/denominator heapや消去fillを含まない。

中心q=0でAの最初42行をC=A_selected(0)とする。保存known-fiveの最初の点は(z,n)=0で、既存certificateのdetC=1680 mod32003は単元である。全QQ係数はp-integralなのでCはQQ可逆、一般betaでもbeta-only minorを反転して同じrecipeを使える。逆行列の新しい展開は行わない。定数可逆行操作

    A'=[ C^-1 A_selected ; A_rest−A_rest(0) C^-1 A_selected ]
       =[ I42+L(q) ; N(q) ]

により、L,Nはhomogeneous linearとなる。b(q)A'=e_jでb=(s,r)、deg b≤1ならs0=e_jであり、各一次係数s_vは

    s_v=−e_j L_v−r0 N_v

から決まる。残る自由なcofactor未知数はr0と45本のr_v、計87×46=4002個である。残る条件はdegree2だけで、最大42×1035=43470本。正規化後の係数heightやfillはここで計算していない。一般degreeDでもrの87*binom(45+D,D)係数を自由未知数とし、sは次数ごとのrecurrenceで決められる。ただしD2でも94047未知数であり、直ちに安いsolverになるとは言わない。

中心での正のmembershipは一般betaへ自動的に持ち上がらない。例えばC(beta)=(1,beta)^T、b=(1,0)^Tはbeta=0で解け、一般点では解けない。Cのrankは一定でcokerもfreeである。実coefficient mapのfull row rank、または全beta係数恒等式を満たすrational solution recipe、あるいは普遍compatibilityを別途証明したsurjectiveな制約空間が必要である。今回のtall mapに単に中心解やcolumn minorがあってもこのgateは満たさない。

既存126×126 cross inverseによる `[z_k I−R_k(n)]` 型のnormal表示とcommutator determinantも先行packetにある。n=0では定数commutatorが可逆で局所的formal inverseを与えるが、全normal上でそのdeterminantを単元にすることは未知残余を除くことになる。新しい多項式left inverseや全fiber長さ5は主張しない。任意次数left inverseが不可能とする新しい理論的obstructionも、今回見つかったものとして主張しない。
