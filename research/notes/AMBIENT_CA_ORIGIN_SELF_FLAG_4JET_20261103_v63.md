# B18原点のCA self-flag四次jet：Gmは二つの重複点へ切られる

2026-11-03。独立監査済み三次jet packetを固定し、新しい四次計算だけを行った。
三次で残ったGmは、四次の必要条件によって **t=±1の二つの二重点**へ切られる。
各点上には先行高次自由度が一本残る。有限四次jetは両符号で存在する。
実際の形式的・代数的pair branch、全次数の非存在、conductorは結論しない。

## 全自由度を保持した実係数計算

CA fixed-fiber C150×20、第一sourceのparabolic断面free10、full relative flag9、
M=E10−2E03、および全150行membership projection πを三次packetと同じ座標で使う。
三次までのgと第一sourceは

    g=I+aM+a²B+a³C3+...,
    P=P0+aP1+a²P2+...

である。三次必要曲線を

    v1=t−1, v0=−7t/4−1+3/(4t), v2=0, t≠0

と**正確なLaurent座標**で表示する。
三次の150×18行列Dはrank15で、三次補正(P2/free10,C3/flag8)に
三自由度(y0,y1,y2)が残る。これらを一つも零に固定しない。

    high18 = −D_pivot_inverse*(三次intrinsic)_pivot + T*(y0,y1,y2)

とした全解を保存した。Tは三次packetにあるDの全kernel18×3である。
四次の新補正(P3/free10,flag-order4/free8)も同じDのimageであるため、
これを消去した全150行のcokerが正確な四次の必要十分な有限jet条件となる。

実5×5行列logarithmは

    L1=M,
    L2=B−M²/2,
    L3=C3−(MB+BM)/2+M³/3,
    L4=−(MC3+C3M+B²)/2+(M²B+MBM+BM²)/3

である。新しい四次flag補正は最後にDで扱うので、ここではL4の自由加法項を零に置く。
M⁴=0である。5×5行列としてexp(log g)=gを次数4まで記号的に確認し、
実150係数のLie表現で同じ有限exponentialを計算した。
従ってLie作用の二乗と行列の二乗を混同した近似ではない。
係数環は独自の疎Laurent多項式QQ[t,t^-1,y0,y1,y2]で、有限jet以外の展開打切りはない。

## 正確な四次ideal

全150行cokerのQQ係数rankは4、独立行は零始まり46,63,72,95であり、その式は

    q1=−y2(t²−3)/(5t),
    q2=3(t−1)²(t+1)²/(2t),
    q3=−(4t²y0+(7t²+3)y1)/(2t),
    q4=−3y2(3t²+1)/(20t).

両idealの包含は次の短い係数証明で正確に成立する。

    y2=(3t/2)q1−(2t/3)q4,
    (t²−1)²=(2t/3)q2,
    4t²y0+(7t²+3)y1=−2t*q3.

逆包含は各qの表示から直ちに従う。従って必要idealは

    J4=(y2,(t²−1)²,4t²y0+(7t²+3)y1)

である。四次条件を満たす先行係数のschemeは

    QQ[t,t^-1,y1]/((t²−1)²),
    y0=−(7t²+3)y1/(4t²), y2=0

であり、二つの二重A1である。
新しい四次補正のkernel三自由度も含めれば、正規化された四次jet全体は二つの二重A4となる。
reduced fiberではy0=−5y1/2だが、**t²−1の二乗を単なるt²−1へ置き換えない**。
微小なtの変位は後続liftingで保持する必要がある。
五変数の小さなQQ GB(q*t−1も含む)も保存したが、idealの証明は上の明示係数である。

## 両符号での有限四次witness

t=±1、y0=y1=y2=0を選び、最後の18補正もDの固定inverseで解く。
t=1では新source P3/free10=(0,0,0,0,0,0,2,0,0,1/3)、flag-order4は零。
t=−1では両者とも零である。
それぞれについて、Lie/logの再利用とは別に、実5×5のgとその形式逆行列を直接作り、
外積gおよび二次monomialの逆座標代入による物理150係数を次数0..4まで計算した。
全係数がLie計算と一致し、全150 membershipが成立することを確認した。
第一・第二sourceの150×5係数と5×5×5群行列をcertificateへ保存した。

## 限定と次の判定条件

この結果はordinary flag separationがa一次で第一sourceがO(a)のchartに限る。
有限四次jetが可解でも、全次数へ持ち上がるとは限らない。
Nの19係数模型やR19のacross-CA/CB conductorを個別CAのconductorに転用しない。
B18原点のCA全odd ideal、実枝、ν方向nilとの関係は未決である。

次はJ4の**非被約構造を保持して**全ての高次自由度を加えた五次の判定、
または両符号witnessからの実際の有限多項式paired-source構成が候補となる。
後続の五次計算はこのpacketには含まれない。

再現:

    source /workspace/cloud-setup/activate.sh
    python code/ambient_CA_origin_self_flag_4jet_20261103_v63.py --output <new-output-directory>

code: `code/ambient_CA_origin_self_flag_4jet_20261103_v63.py`。
reportと再利用可能な全150係数:
`reports/AMBIENT_CA_ORIGIN_SELF_FLAG_4JET_20261103_v63/`。
制限は120秒・2.2GB、実行は1秒以内。大規模GBは行っていない。
