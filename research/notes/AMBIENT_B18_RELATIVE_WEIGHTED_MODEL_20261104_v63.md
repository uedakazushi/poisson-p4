# Actual B18 full D を零weight17変数上の三変数三次以内へ戻す

2026-11-04。P0=e_y∧(x²∂h+hx∂z+(h²+z²)∂w)と、固定torus diag(1,1,λ⁻¹,1,1)・tensorλを使う。対象は元の非被約D=CA∩CB。前packetの全20方向positiveという方法は使えないが、相対gradingは残る。

**結論**：共通G19横断完成環は、B=C[[X1,…,X17]]に対して

    D_perp ≅ ( B[Y1,Y2,Y3]/K )^(m_B,Y)-completion

と書ける。KはY次数に関してhomogeneousで、次数0,1,2,3の式だけで生成できる。各個別成分のminimal20変数kernelは次数≤2だが、全Dでは元ambientのgraph差を保持するため **≤3** が安全な上限である。可能な三次graph差は6本。各次数のmonomialは1,3,6,10個なので、実係数の計算対象は計20種類のmonomialに付くformal零weight17変数函数へ限られる。これらの係数は未計算である。

## 1. 全次数でのsourceと共通切断

既受理の二つのnormalizerは、同じ29個のphysical target函数R29上でfinite free rank2である。保存された源core10の全150 polynomialを再読し、p,k0,…,k8中でk2だけをweight1、残りをweight0としたとき、各出力が対応するphysical weightにhomogeneousであることを全項確認する。これは一次jetの推測ではない。源kernelのLie matrixはtorusと可換で、源速度のweight1 fiber成分も零である。従ってkernel pのcharacterは0である。

source core10と19群方向のnormalizerへの写像は、既存の29次unit blockから形式的etaleである。群sectionはscalarとGLのweight eigenvectorsで選べ、torusの共役で保たれる。sourceの上の明示的weight作用と全polynomial source写像はcompatibleである。formal inverseの一意性により、この積分解も全次数でtorus-equivariantとなる。源regular normalizerの全characterは0^24+1^5、群19を除くと0^9+1^1である。有限・birational normalizationの群作用への一意liftとも一致する。

共通R29の中から、保存されたorbit19 matrixのunit minorに対応する19個のphysical函数を切断に選ぶ。そのweightは0^15+1^4、残る共通R10は0^9+1^1である。双方の成分・normalizer・conductor・B/A moduleに同じ作用積を適用する。したがってこれらのmodule上で19切断はformal変数の正則列となる。特にA→normalizerの単射は切断後も保たれる。単にfinite mapを任意にbase changeして単射を仮定する議論ではない。

切断後もnormalizerは同じR10上free2かつregular10である。唯一のpositive base parameterをw、零weight9変数環をR0と書けばR10=R0[[w]]。trace-zero rank1 moduleのfiber characterはpから0であり、linearly reductive torusのequivariant Nakayamaでweight0のgenerator τを選べる。従って

    B_i=R10⊕R10 τ_i, τ_i²=ν_i∈R0,
    A_i=R10⊕I_i τ_i。

ここでI_iは**actual common base上のideal**であり、新しいabstract99モデルのbaseをCA/CB間で同一視していない。保存された全actual odd cotangent basisによりI_i/mI_iのcharactersは0^8+1^2である。equivariant Nakayamaで八本f_j∈R0と二本w*g_k (g_k∈R0)を選べる。よってI_i=I_i0 R10+w J_i0 R10となる。regular sourceのtorus-linearizationは明示した全source mapとformal etale座標に基づき、微分だけから高次を捨てていない。

## 2. 個別成分のrelation bound2

I_iの八本degree0・二本degree1生成元のhomogeneous syzygyを考える。degree0は八本f_jのR0上のsyzygyである。degree d≥1の任意のsyzygyは、八成分がw^d*c_j、二成分がw^(d−1)*d_kという形なので、regular domain R0[[w]]内でw^(d−1)を除けばdegree1のsyzygyになる。したがって全first syzygy moduleはpositive degree0または1で生成できる。R0上の二つの有限syzygy moduleから有限生成も従う。

odd target Z_j=b_j τ_iを使えば全algebra kernelは、このfirst syzygyと

    Z_j Z_k−ν_i b_j b_k

で生成される。後者のdegreeは0,1,2である。pair productsで任意の式をR10+sum R10 Z_jへ還元した後、basis(1,τ_i)で残る関係がちょうどfirst syzygyになるため、全kernelに対する主張である。

元のphysical target座標に戻るeven shiftsも同じweightを持つ。共通minimal target20は零17・正3なので、weightを保つambient formal automorphismはX→φ(X)、Y→M(X)Yという形に限られる。したがって個別minimal kernelのdegree≤2は保たれる。既受理99関係によるgenerator上界もminimal Cohen presentationのformal一意性とequivariant Nakayamaで移る。各個別kernelは高々99本でよい。abstract99をw=1,その他base=0,Y_wr/Y_ws=1でweight付けした64/32/3の分布は一致検査に使うが、その新しい座標をrelative位置と同一視しない。

## 3. 全ambient graphを残すのでfull D boundは3

canonical ambient126のcharactersは0^36+1^60+2^24+3^6。共通actual tangent39は0^32+1^7である。従ってminimal39座標以外には、weight0/1/2/3それぞれ4/53/24/6個のgraph座標がある。全群切断後もこの87個は残る。

それぞれの成分について、linear kernelの87 unit blockをequivariant formal IFTで解けば、graph座標は共通minimal座標のhomogeneous函数Φ_iである。零weight部分はformal零変数函数、weight1はpositive変数の一次式、weight2は二次、weight3は三次となる。有限jetで打ち切ったからではなく、positive weightとequivarianceから**そのweightを越えるmonomial自体がない**ためである。

CAのgraphを使って87変数を消したfull D idealは、個別二kernelに加え全87本のΦ_CA−Φ_CBを含む。別々の99モデルをrelativeに揃えたとは仮定しない。このgraph差を落とすと元のnilも失い得る。全式はdegree≤3で、上界本数は99+99+87=285。degree3が必要となり得るのはambientの六個のy²∂a∧∂b(a,b≠y)だけである。これら六差のlower-degree idealへの所属は未検査であり、whole D degree≤2は主張しない。

## 4. 零weight fixed schemeとnil class

非負gradingなのでpositive3変数を零にしたquotientはdegree0 subringそのものである。

    D_fixed= C[[X17]]/K0,
    K0=J_CA,0+J_CB,0+(四本のzero-weight graph差)。

個別成分のfixed ringはexactにR0⊕I_i0 τ_iであり、normalizerは零weight9変数のregular ringである。actual係数I_i0や二つのrelative像は未計算である。

前に独立受理されたactual Dのnil-cotangent basisはu63(weight0)とu46+u95(weight1)である。nil idealはtorus-stableかつ完備閉なので、u63に対応するnil liftのweight0成分もnil idealに入り、同じ非零linear classを保つ。degree0部分はDのsubringなので新しいnilは作らず、(m_D²)_0=m_fixed²でもある。従ってD_fixedのedimは17、そのnilradicalのlinear cotangent像は**正確に一次元**でu63が基底である。生のu63自身のnilpotenceを主張しない。

positive3がregular sequenceだとは分かっていない。このfixed quotientのdepthやnilをfull Dのdepthへ移さず、full nil ideal・socle・depth19/20は未決のまま残す。

このpacketは全source polynomialのweight/sparsityと既存finite-source kernelを結ぶ有限の構造定理である。係数の消去・GB・新点・primary solverは実行しない。
