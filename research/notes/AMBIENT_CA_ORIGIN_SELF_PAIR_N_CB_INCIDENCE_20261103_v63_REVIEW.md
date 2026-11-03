# CA self-pairの実N/CB incidenceと一般点でのNの滑らかさ

2026-11-03。個別reduced CAの二旗と、個別reduced Nの全旗fiberを区別する。

## 全次数の物理source監査

raw CAを独立にcanonical化して得た150×20行列Cから第2列（0起算）を除くと、保存されたNの150×19行列と全係数で一致した。
同じ旗のparabolic16個のLie作用を実5×5テンソルで独立に計算し、全304列がN19に入ることも確認した。
連結parabolicの不変部分空間なので、既知の全N proper模型は全CA proper模型の**同じ旗base上のclosed rank19 subbundle**である。
これはS61だけのsourceではなく全射影模型についての包含である。

第一表示P(a)のC[2]係数は恒等的に0、第二表示g_*P(a)では−a²である。
従って全aでP(a)は実N sourceに入り、a≠0で第二CA旗はN sourceに入らない。
実Nの19係数も全150行で再構成した。

またGの第一行を(−3a,−1,0,0,0)、他行を恒等とするとdet G=−3aであり、

    f_B=3a x²,
    H_B=−a xwz+xw(h+x)/(3a)+(h+x)²z/(9a²)+z³/3

を持つraw CBのcanonical代表とG_*P(a)がQQ(a)上の全150係数で一致する。
これは全次数の有限座標変換であり、標本やjetからの外挿ではない。
元のPと変換後のテンソル双方について25 divergence係数と350 Jacobi係数も独立に零と確認した。
a=0でのCB membershipは既知の共通P0、または閉像の閉性を用いる。

## 個別Nの全fiberとsmoothness

別の独立監査 `AMBIENT_CA_ORIGIN_SELF_PAIR_GEOMETRY_20261103_v63_REVIEW` により、
a=0を含むtarget openへ全CA proper射を制限するとfiniteかつ全fiber length≤2であり、
その逆像の一般のa≠0では全CA fiberが正確に二つの被約点となる。
全N sourceは上のclosed subbundleなので、その全fiberはこの全CA fiberのclosed subschemeである。
一方の旗はNに入り、第二の旗は−a²≠0によりNに入らない。
従ってそのopen arcでは**全N fiberが正確に一つの被約点**である。
他の旗を同じ局所flag chartの外で見落とす推論はしていない。

N19とflag9から独立に作った実150×28微分はQQ(a)でrank28。
これはCAの全29方向から一つのfixed-fiber方向を除いたものでもあり、a≠0では単射である。
N sourceはsmooth28、source→targetは上のopenでfinite。
閉像を定めるreduced N局所環Aは有限なsource局所環Bへ単射であり、全fiber length1からB/m_A B=C。
従ってB=A+m_A Bであり、有限A加群B/AへのNakayamaからB=Aとなる。
よって一般のこのarc上で個別reduced Nは**cone smooth28、射影smooth27**である。

これは全a≠0でのfiber定理、CBの滑らかさ、CA原点の全odd ideal、全conductor、actual CA/CB intersectionの型を主張しない。
二つのCA旗から二つのN旗を推論せず、既知のN subintegral completionとの矛盾もない。
R19のmembership・非membershipはここでは分類しない。

## 再現

    PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/ambient_CA_self_pair_N_CB_review_20261103_v63.py --output reports/AMBIENT_CA_ORIGIN_SELF_PAIR_N_CB_INCIDENCE_20261103_v63_REVIEW.json --replay work/v62-review/CA-self-pair-N-CB-incidence-replay

著者のfresh certificate byte一致、全physical source恒等式、parabolic不変性、Nの正確なrankとminorを保存する。
全proper模型の被覆と全CA fiberのopenは既存の独立監査を再使用し、GBを再実行しない。
