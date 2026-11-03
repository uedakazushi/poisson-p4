# 実全四次の122×64判定と90×64十分条件

2026-11-03。既存の全四次319×239と286×238の十分条件に、heavy-linear elimination と pure-quartic quotient を加える。**元の全33本の genuine quartic を保持する122×64判定** が新しい点である。係数の大きな展開・GB・有限体・新点検査は行わず、実在sourceを指定した有限 circuit と可逆操作を保存する。rank64、全四次全射、長65は未証明である。

基礎開は U0={aFBC!=0}。追加の chart 条件は、16 genuine cubics の heavy7 projection の7minor `DeltaH!=0` と、168×210 pure quadratic-products の146minor `DeltaQ!=0`。これらは追加 chart 単元であり、全 core の既知 unit と呼ばない。FULLU={adFBCJ!=0}では既証 product41と、今回別に証明した quadratic-product146の有限欠損結果に接続するが、本縮約は full cubic57を仮定しない。

## 1. 実際の source

七 light を x0,...,x6、heavy v5 を H と書く。xの元座標順は `[v0,v1,v2,v3,v4,v6,v7]`、重みは1と2である。

既存の厳密自由57行基底は、41 product行と16 genuine cubic行 G_j からなる。57×91の全係数は既存の圧縮sourceに固定済みである。この自由性は周囲91空間への点別単射性を意味しないが、全三次行の生成元としての等号は U0上で保たれる。

全重み四の同じ実像は、既存319 packetにより

`6 quadrics *29 +16 genuine cubics *7 +33 genuine quartics`

で生成される。元135式のweight4非選択33行を一つも捨てない。circuit は全33 source indices と weighted implicit graph oracleを固定する。元135式・選択34行・厳密normal chartは先行packetと同じで、新しいquartic sourceを作らない。

## 2. Heavy-linear cubicsによる28列の消去

DeltaH chartで七つの G_E を選び、その H*x_i 係数行列を T とする。`f_i=sum_k(T^-1)_ik*G_Ek` と取れば

`f_i=H*x_i-p_i(x)`

である。p_iはpure cubic。残る九行は heavy coefficientを引いて

`r_j=G_j-sum_i coeff_(Hx_i)(G_j)*f_i`

とし、pure cubicになる。16行のこの変換はDeltaHだけを逆にする可逆変換である。

49行 f_i*x_j の H*Sym²light 係数は、各行がただ一つの unordered pair {i,j} に1を持つ49×28定数行列である。i<=j の28行が恒等行列となる。残る21行との差は

`f_i*x_j-f_j*x_i=p_j*x_i-p_i*x_j (i<j)`

というpure quartic。従って全49行は28個の単位pivotと21 commutatorsへ正確に置き換わる。codeは定数mapのrank28、21本のkernel basis、全indexを検査する。

以後のnormal formは `H*x_i*x_j -> p_i*x_j (i<=j)` と固定する。六行H*q_iはこの規則でpure quarticになり、九行r_jの七light倍は63本のpure quarticである。この段階の低重み residual は正確に `6+21+63=90` 行。

## 3. Heavy-squareは全開で一列消去できる

実genuine quartic Q77,Q131の H²係数を b77,b131とする。既存の exact QQ certificate は u,v に対して

`u*b77+v*b131=U=(50aBC)^2`

を与える。個々の b77,b131を全開で逆にしてはいけない。その代わり二行を

`chi=u*Q77+v*Q131`,
`psi=-b131*Q77+b77*Q131`

へ変える。2×2変換行列 `[[u,v],[-b131,b77]]` のdeterminantはUという既知単元。psiはH²係数0、chiはH²係数Uである。他の31 genuine rows Qkから `(b_k/U)*chi` を引けばH²を消せる。先のH*Sym² normal formの後、chiから

`H² -> -NF_H(chi-U*H²)/U`

を得る。こうして33本は一つのH² pivotと**32本のpure quartics**に正確に置き換わる。元の33本全てのspanを保つ。追加の b77!=0/b131!=0 chart は不要である。

codeは保存されたactual b77,b131,u,vのQQ恒等式、普遍2×2変換を照合する。別の小fixtureでは未知hがpsiに残ることも確かめ、未知因子の除去を認めない。

## 4. Pure quadratic quarticsの146次元を割る

六quadrics×28lightquadraticsは168行、pure light degree4は210列。DeltaQ chartで146行I、146列Kを選び、P=Q[I,K]を逆にする。補集合Jは64列。pure quartic row zのquotient coordinateを

`pi(z)=z_J-z_K*P^-1*Q[I,J]`

と固定する。

先行の全点upper146により、残る22 quadratic-product行のpi像は零。これはrank147 minorが恒等的に零というalgebraic upperを用いるため、単なる標本点での消滅ではない。DeltaQ chartで146方向の可逆消去を実現する。

したがって元239列には28 H*Sym²、1 H²、146 pure quadraticの合計175単位pivotがある。残る64列への実行列の行は

- 6行の NF_H(H*q_i)、
- 21行の commutators、
- 63行の r_j*x_i、
- 32行の genuine quartic residual。

全てにpiを適用する。この **122×64** 行列を C122 とすると、任意のgeometric pointで厳密に

`rank(actual full weight4)=175+rank(C122)`。

source generation・可逆row変換・Schur操作から、単なるrank等式に加えてcokernelもこの64方向のcokernelと同型である。従って全四次239方向の全射は `rank(C122)=64` と同値である。

## 5. 90×64は十分条件

最後の32 genuine quarticsを省いた C90 は最初の90行である。低重みproductsはH²係数0だから

`rank(low-weight products)=174+rank(C90)`。

rank(C90)=64なら低重み238方向が全射となり、global H² unitを加えて全239方向が全射となる。これは十分条件であり、C90のrank不足を元のfull quarticの不足に読み替えない。32本のgenuine residualが不足を補う可能性を保持する。

これで既存286×238十分条件と同じ低重みimageを、追加chart上で90×64へ表した。新しい122×64は全319×239と同じimage/cokernelを扱う。full cubic57はどちらの縮約にも不要である。

## 範囲・次の有限目標

DeltaHとDeltaQの全minorを取れば、heavy projection7かつpure quadratic-product146の開を有限に覆う。FULLUの任意のcore曲線の一般点では、quadratic-product146は今回の有限欠損定理で保証される。heavy欠損曲線は既存の別proper certificateで扱われており、そのscopeを任意の曲線へ外挿しない。

現在未証明なのは、これら実122×64行列のrank64不足に曲線があるかどうかである。90行の十分条件だけを試す場合も、失敗を元full imageの失敗とはしない。全三次57、全四次239、横断長65、ambient除去は依然として別の主張である。

再現は `python code/s61_nonsquare_quartic_122_circuit_20261103_v63.py --output <new-output>`。元圧縮57×91sourceのhashと全112heavy成分、33quartic source labels、定数49×28/kernel21、actual H²unit、全chartとprojection定義を検査・保存する。大きなparameter係数は展開しない。
