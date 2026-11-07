# 実際の最大添字 H と raw 射影による Kbar 二乗の消滅検証

独立 SOURCE 検証は PASS。実際の original-z ordered prefix は最大の有効添字を選ぶ。これを使って全ベクトル H を再構成し、sector の指数 (1,1,2) と、最後の sector の `w_i=v_i+r²H_i` も残してから次数二を抽出した。U、S の一次係数や八つの Qeven 一次係数は代入していない。

36 二次単項式と九つの形式方向について、最終三次元射影も、実際の三本の二次 relation を加える前の既知174次元部分空間の商六次元射影も、全324列が零である。実際の積に必要なのは、L_U、L_S の r 係数が零であるため、zz の28単項式に対応する252列である。残る72列では、cycle でない形式入力上の pure-r 最終成分を零へ延長する規約を採っただけであり、任意の非 cycle に対する contraction を主張しない。必要な252列では、この延長の選択は使わない。

凍結済み leading packet の min-active という中間記述は、実際の shared H の字義的な実装とは異なる。その記述を修正したと仮定せず、ここで実際の max-active 演算を別コードで検証した。元ファイルは変更していない。

raw128 入力 `(free44=0,H2)` の Pi1 も独立に導いた。直接 section は `n=E(mE)^−1` で E の42行が零なので、n の42行は全次数で零である。最初の24 unit 消去および次の18 unit 消去はこの入力をそのままにする。元780列の X1 から81個の residual first を再構成し、元 W0 minor の保存 inverse の左右両積を全成分照合すると、raw Pi1(head0,r) は元 H27 成分に **+8/3** を掛ける。RIGHT106 の relation 側の −8/3 を raw 入力へ移用していない。旧 edeef の符号誤同定は保留のままである。

実際の三次 raw 積は `L_U L_V(q Wq+B Wbeta)` と e43 の交差項からなる。m0e43=0 であり、m0 の sector0 最終行は Wq、Wbeta の両方を零へ送る。従って whole H の H27 次数二は恒等的に零で、正しい raw 係数 +8/3 の寄与も零である。

既知174部分を消す際には、一次の斉次 relation だけでなく真の全関数35 relation を用いる。この真の M2 の既存四方向 callback と fixed174 section の r 係数を独立に合成した結果は、head0 の z0r、z4r に係数50/3、125/3をもつ二項だけである。raw 二次入力の free44 は零、Pi0 第一行は e42 なので、その両係数は零である。従って cubic clearing を捨てているのではなく、真の補正 functional 全体が零を与える。

L_U、L_S に r 成分がないため raw42 の pure-r 三次係数も零となる。高次 Gamma、section、m1、Pi の未知項を零と置かずに次数を追うと、この pure-r 三次 head0 dual に入るのは上述の三項だけである。真の relation で二次部分を消した積に r を掛けると次数四から始まり、完全1650座標の四次 dual は head0 r4 の一座標だけを読む。m5Tnil=0 により、実際の `rU²=rUS=rS²=0` が従う。

受理済みの q2 による bar 消滅と全 r-closure から、全体について `r Kbar²=0`、`Kbar²⊂m³Tnil` となる。既知の m²Tnil の r-Jordan 型は (3,2,1) で、m³、m⁴の長さは3、1である。従って r m²=m³、r m³=m⁴であり、m³上の r kernel は二次元。よって **dim_QQ Kbar²≤2** である。

この結論は bar 上の積に関するものだけである。D-finite torsion saturation により H0/ellH0 は K/ellK に単射し、N/K が D-free normalizer nil に入るため K/ellK は Nbar に単射する。ell²K²=0 より K²⊂H0 であるが、生成元数が少なくとも三の H0 を、二次元以下の Kbar² の像から D 上生成することはできない。H0 の長さ、ell 消滅指数、元の nil 三乗や nil 指数、積三本の実際の全値または正確な階数は求めていない。

独立コードは producer の H・行列演算を import せず、標準 Fraction と既存の独立 JSON field reader だけを用いる。SOURCE 数値所有上界65,309、入力1,200,071 bytes。事前高さ上界は tensor257、true clearing127、raw C1 313 bits、実測 stored89／pre156 bits。60秒、2GB、150,000数値、stored384／temporary1153 の範囲内で完了した。新しい積値、物理係数、jet、点、素数を要求していない。
