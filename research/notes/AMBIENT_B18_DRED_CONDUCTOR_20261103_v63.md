# generic B18 の元の被約交叉：実際の生成元による全 conductor

2026-11-03。対象はactual D_redの形式的完成環であり、元の非被約Dのconductorやdepthではない。
一つの承認済みN/R二変数四次source計算だけから、欠けていたT³p係数を求める。
全28変数係数の巨大な展開はせず、実際のsourceで一意に定義される生成元と、有限な有理数係数で
**全conductor ideal**を指定する。抽象normal formとの未確認の座標同一視は用いない。

## 1. 元の実際の生成元を定義する全source回路

既受理の物理N19固定空間、GL19方向、補償されたkernel方向Kを用いる。
全sourceは `G_outer=I+Σg_i E_i` と `G_inner=exp(pK)` の作用である。
固定sourceは `P0+Σu_j C_j−p L_K(P0)`、C_jはC2を除いた固定8方向である。
これらの全150係数、GL添字、K、補償列をsource証明書に保存する。
outerの逆は原点で単位の有理行列であり、この式は任意次数の実際の全sourceを指定する。

27個の実際の標的座標をregular base変数として固定し、sourceの27変数を形式的陰関数で一意に解く。
残る実際のfold行27を28番目のbase座標とする。
fold方程式の単位の定数は−3で、distinguished monic二次式を

    p²−c1(R)p−c0(R),  tau=p−c1(R)/2,  nu=c0(R)+c1(R)²/4

とする。従って正規化は実際の物理base R28上でfree rank2、`B_N=R[tau]/(tau²−nu)`である。
全physical出力をこの二次式で割って `q_j=a_j(R)+b_j(R)tau` とする。
選択されたactual最小N37座標の後ろ九行

    [18,19,28,29,60,61,67,46,72]

のb_iを用いる。この九つは元の全奇係数ideal Iを生成する。
実際、選択37座標は全完成像環を生成し、二次cover内でodd積もこれらb_iのR倍に戻るためである。
よって `A_N=R⊕I tau` は元のactual像であり、normalizer内のconductorは `I B_N`。
そのR-module生成元は **実際の**18個 `b_i,b_i tau` で尽くされる。
本稿のb_iは抽象的なs_iではなく、上の全source/IFT/Weierstrassにより一意に定義される関数である。

## 2. 承認された一回のN/R source4

同じ保存済みの実際のR曲線を標的とし、全次数4の(T,p)だけを用いる。
各次数1–4の27IFT補正と左右逆行列、全150source係数を保存する。
計算は有理数添字Lie演算で `log(I+M)` とkernel指数を四次まで作用させる。
元の25divergence、350Jacobiを二変数全次数4まで直接検査する。
全次数3以下の係数は旧独立受理済み証明書と完全一致する。

fold方程式を単位×monic二次式として四次まで掛け戻す。
全150出力を割り、b_j(T)をT³まで得る。これが以前欠けていたT³p係数である。
tauのtrace shiftも行い、元のC4=k[T]/(T⁴)への実際のgluing写像phiについて

    sigma_i=phi(b_i),
    eta_i=phi(b_i*tau)=q_Ri−phi(a_i),
    phi(nu)=64T²

をC4の全係数で求める。ここでtauそのもののphi値を仮定していない。
phiは既受理の正確な接触k=4から存在するA_N→C4である。
normalizerのfree rank2座標をbase-changeして係数を評価するのであり、非平坦な像idealのkernelを
base-changeで交換する議論ではない。

## 3. 元の被約交叉のconductorの二条件

同時G19因子を取り除いた被約環は `F=A_N×_C4 k[[T]]` である。
正規化は `B_N⊕k[[T]]` となる。
二つの正規化成分を別々に掛ける条件から、そのconductorは正確に

    J ⊕ T⁴k[[T]],
    J={a∈I B_N : phi(a)=phi(tau*a)=0}.

B_NはA_N-moduleとして1,tauで生成されるので、この二条件だけでaB_N⊂ker(phi)が尽くされる。
tau*aはa∈I B_NならA_N内にあるため、二番目のphiも定義されている。

従ってPsi(a)=(phi(a),phi(tau*a))とすると、その18列は

    b_i       ↦ (sigma_i, eta_i),
    b_i*tau   ↦ (eta_i, phi(nu)*sigma_i).

Tの0–3乗倍を含めた **8×72のQQ行列** を保存し、その全核68次元を有理数で計算する。
各核ベクトルの恒等式と自由列上の独立性を検査する。
このQQ核だけをidealと呼ぶのではなく、次のbase-curve kernelも加えて全moduleを復元する。

## 4. 一つの曲線とG19を混同しない

全coneのN正規化は28次元だが、Rの全cone成分は20次元である。
一つの物理R曲線だけを全R成分と呼ばない。
偶base上の19軌道方向ORから、非零な横断R曲線速度を捉える線形式ellを選ぶ。
`ell*OR=0, ell(v_R)=1`を有理数で検査する。
保存された**正確な**物理R曲線をr(T)とし、

    u = [ell(r(T)) の形式的合成逆関数](ell(R))

と定義する。u(r(T))=Tであり、duは19軌道方向で零である。
19個の軌道pivotの偶base座標について

    h_j=R_j−r_j(u)

を取る。これは切断されていない実際の形式関数である。
h_jの軌道微分は元の単位minorそのままで、全R曲線上で零となる。
従ってh19=0はこの曲線を含む実際のtransverse sliceであり、受理済みの群作用の形式的IFTを適用できる。
N、正規化、conductorの対、R、被約交叉を同時にG19×sliceへ分解する。
偶baseで切るため、二次coverのbasis1,tauも保たれる。

base-curve写像rho:R28→C4の核も有限式で指定できる。
`ell(r(T))=T+bT²+cT³+O(T⁴)`なら
`u3=ell(R)−b ell(R)²+(2b²−c)ell(R)³`とする。
ellの非零係数の座標を一つ省き、残り27個の
`R_i−r_i(u3) mod u3⁴`とu3⁴でker(rho)を生成する。
この核と、68個のQQ核ベクトルのT^kをu3^kで持ち上げたものから、

    ker(Psi)=ker(rho)*I B_N + Σ_R(持ち上げたQQ核)

という全module kernelが得られる。
19個の真のh_jでsliceを取り、その後G19を戻したものが元のconeのconductorである。
全R28の一曲線核そのものをcone conductorと同一視していない。

## 5. cyclic image、完全交叉長さ4、17生成表示

物理行61の偶生成元b_mを取る。保存された係数では

    sigma_m=46T+O(T³), eta_m=232T²+O(T³).

g=Psi(b_m)と置く。tauの作用は `(a,b)↦(b,64T²a)`である。
`g,Tg,tau*g,T²g` の8×4行列に単位4次minorがあり、全Psi像を張る。
従って `I B_N/J` はB_N-moduleとしてb_m一つで生成され、長さ4である。
九個全てについて保存した有理数係数により

    Psi(b_i)=h_i(u,tau)g,  h_i∈span_QQ(1,u,tau,u²)

となることを全8成分で検査する。

mu=116/23と置く。直接の有限係数検査により

    (u*tau−mu*u²)g=0, (tau²−64u²)g=0.

`mu²−64=−20400/529≠0`なので、この二つの二次式からu³が従う。
実際、`u*q2−tau*q1−mu*u*q1=(mu²−64)u³`である。
残る偶base座標は、実際のnu,uと26個の物理座標で正則座標を作る。
この全28次Jacobianの単位minorを検査する。
26個の座標c_jについて `c_j−r_j,1*u−r_j,2*u²` が像を消す。
T³は像を消すため、ここに三次以降の曲線係数は不要である。

これら26個と二次式二つのidealをK_auxとすると、正規化の商は正確に長さ4で、
像の長さも4であるから `K_aux=Ann_BN(g)`。
19個の正確なgroup gaugeはこのidealに入り、sliceでは七つの偶一次制約と二つの二次式、
計9生成の完全交叉ideal Kを得る。
そのArtin環のHilbert functionは `(1,2,1)`である。
一方、R-moduleとしてのPsi像はT作用の階数 `(4,2,1,0)`を持ち、Hilbert functionは `(2,1,1)`。
この二つのfiltrationを混同しない。

全J idealは明示的に

    J=( b_i−h_i(u,tau)*b_m (i≠m),  K*b_m )

で生成される。8個の差とKの9生成元の積、計17個の**生成上界**である。
証明は、任意のa∈I B_Nをこの差で `s*b_m`へ戻し、Psi(a)=0がs∈Kと同値であることによる。
これは `J=I∩K` という別のideal同一視ではない。最小生成元数17とも主張しない。

## 6. 保存結果の範囲

元のD_redのnormalization conductorは、上記actual-generator JとR側T⁴の組で全体が決まり、
G19を戻す。N側の元のnormalization conductor I B_Nとの差 `I B_N/J` はsliceで長さ4、
G19復元後は19変数の正則係数環上rank4の自由moduleとなる。
全conductor quotientは有限長環ではない。実際 `F/conductor≅A_N/J` であり、
既知のN conductor quotientへ有限長4の部分moduleを付け加えたものとなる。
その環の直和分裂や、全28変数の展開済み多項式表示は主張しない。

source4と全18写像の実係数は `source4.json.gz`、8×72核、cyclic係数、実際のgauge、
全Jの17生成recipeは `conductor-kernel.json.gz` に保存する。
これらの生成元は上の全sourceの一意な形式回路に束縛されているため、有限jetだけのconductor候補ではない。
新しいsource次数はこの一回の全二変数次数4のみ。60秒・2.2GBを超える場合は停止し、自動拡大retryしない。

* code: `code/ambient_B18_Dred_conductor_20261103_v63.py`
* report: `reports/AMBIENT_B18_DRED_CONDUCTOR_20261103_v63/`
