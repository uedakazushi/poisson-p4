# q2 の genuine Nbar 消去に対する独立 SOURCE 査読

primary `ee0ec8dc…` と specific H companion `eb63cdc8…` を査読した。独立報告は `reports/AMBIENT_B18_C_KBAR_Q2_BAR_ANNIHILATOR_SOURCE_20261105_v63_REVIEW/report.json`、SHA-256 `caed49bc0049782a568e30b738068869d864395d8ef754be6edfdb4d69606f6d`。結論は **q2·Nbar=0、従って g·Nbar=0 の bar 上での証明 PASS**。ここで Nbar=N/ℓN、q2=i(r² oldh4)、g=rq2 である。元の R0 上の lift による消去は主張しない。

元の product44 が用いる Vrad と、honest モデルの Vgamma=P⁻¹Vrad を全3×44参照で同じ実関数に結び付けた。m は同じ regular S 行列で、S=QQ[[zA0,…,zA6,rA]]、d1m=Vgamma である。canonical section は minimal node8 の E(mE)⁻¹であり、以前の npre を代用していない。点 m0n0=I と E の支持も独立に検査した。

自然写像の実際の origin 列と、minimal→free44 の定数 section を合成すると、oldh4 の free44 origin は e43 となる。Rbar 線形性から、q2 の選んだ代表は r²e43+O(m³)。したがって a36 と bz7 は m³、bgA は m²に入る。一方、genuine Nbar の実際の5生成元は、最初の origin が e43、残りが零であるため、任意の Nbar class に a36,bz7∈m、bgA∈S を満たす genuine cycle 代表を選べる。これは任意の Tnil class に関する主張ではない。

全36 Qrad-even は proper J⊂m に属する。従って nil 代表 n の b0 は m²、q2 の b0 は m⁴。CA の最初7関数が m、真の gA505938 が全混合方向で m³に入ることから、q2 の B=Σfi bi は m⁴に入る。full product44 の実 DAG を min-plus の次数半環で独立に辿り、全44出力を確認した。結果は

`a_product∈m⁵, bz_product∈m⁵, bgA_product∈m⁴`。

trace monic の全乗算は regular であり、次数を下げない。したがって raw free44 の積 c の次数4成分は e43 だけにあり、c42∈m⁵。実 m0e43=0 なので mc∈m⁵となる。実際の積係数は一つも計算していない。

ここで使う H は、honest の literal 共通 template による三 sector `(z7,r),(z7,r),(z7,r²)` の一回の操作である。全24本の d1 と84本の d2について、符号・power・順序を元 resolution と照合した。任意の真の formal cycle v∈m^d に対して、z=0 の cycle 恒等式と r^n の非零因子性から v7(0,r)=0。ordered Hadamard による Hi は zi を一回割るだけなので m^(d−1) に入り、wi=vi+r^nHi は n≥1 より m^d に入る。

次の pairij は `(Euler_z+2)⁻¹(∂i wj−∂j wi)` である。一般 monomial z^αr^s e_j に対し、δz dz+dz δz=(|α|+1) を49成分の α に関する形式一次多項式恒等式として独立に確認した。分子の z 次数は |α|−1 なので、template の Euler+2 の分母は |α|+1となる。この係数ごとの除算は普通の m 次数を変えず、z 微分一回だけが次数を一つ下げる。よって **全 formal cycle に対し H(m^d cycle)⊂m^(d−1)**。r や r²で割る操作ではなく、原始 jet を得るための query-loss とも別である。

この H を積の**全関数 mc**へ、有限 NF より前に適用する。780列で選んだ H を S 線形に延長して任意の積へ適用したり、係数を H の外へ出したりしない。mc∈m⁵なので Hmc∈m⁴。raw128 の iota(c)=(c,Hmc) は m⁴に入り、d2Hmc=mc から、どちらの right section n についても

`Phi_n iota(c)=c−nmc+n d2Hmc=c`

が厳密に成り立つ。全44行の形式相殺を確認した。section transport J_n(s,b)=(Phi_n(s,b),b) もこの iota を固定するため、raw128 の free44/H84 の領域が変わって e42 補正が隠れることはない。

完全 minimal 射影の Pi0 第一行は `[e42,0_84]` で、H84成分を全く読まない。c42∈m⁵、iota∈m⁴なので、Pi0 の第一成分の次数4は零、regular Pi1以降も次数5からである。完全1650成分の quartic dual は最小生成元0の r⁴を読む entry329だけだった。よって積の m⁴/m⁵ class は零。受理済み Tnil の m⁵=0 により積の Tnil class は零となる。Nbar→Tnil の genuine な単射と、実 D0 積を運ぶ product44 により、q2·Nbar=0 が従う。Tnil 全体に積を与えたわけではない。

受理済み Kbar=(U,S,q2)_Rbar と組み合わせると、Kbar² は Rbar 加群として U²,US,S²で生成される。QQ の span や作用の閉包を調べるときには、特に r 倍を保持しなければならない。三つの積の独立性・値・閉包次元は計算していない。

全参照ハッシュ、実フレーム、H template、matrix domain を検査し、owned 上界107,779、実行0.524秒で完了した。新しい source、IFT、gauge、jet、符号付き積は求めていない。primary と既存報告のバイト列は変更していない。元の q2 lift の上階消去、ℓの取り消し、H0、nil index の改善は結論しない。
