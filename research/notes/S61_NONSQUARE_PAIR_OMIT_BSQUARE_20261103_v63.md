# B²列を除いても actual pair は全域 rank2：quad78 への小さい gate

2026-11-03。actual pair の空間は `A=X=W=0` であり、B は自由である。先行全域QQ pair theoremと同じ自然core座標を用いる。2×6係数列の順序は `(B²,B*vlight,vlight²,B*u,u*vlight,u²)`（既知 unit によるB座標スケールを含む）。今回、**B²列だけを除いた2×5も全 a*P*J!=0 上で rank2** であることを証明した。

## 1. 先行 finite candidate の literal QQ cleanup

f=(v+3)(3v+4)!=0 の場合は、先行 proof が用いた三つの必要条件（列2,3,4とu²列のminor）がすべて今回も残っている。従って C/Q23 の場合分け、C枝の除外、Q23枝での Hbar(v),tbar(v) の決定は、先行の全恒等式をそのまま再検証して使える。

Q23上では、今回残る **第1列とu²列の actual minor** を評価した。これは保存多項式 m1(v) となり、保存された別多項式により `m1*inv1=1 mod Q23` が厳密QQで成立する。従ってB²列を使わずこの枝を除外できる。

f=0 の二枝でも、第1行の残した列は零にならない。実際、`row0_col3-(4/15)row0_col4` は v=-3 で `3t/625`、v=-4/3 で `41t/20000` という既知 unit となる。rank<2ならQ3のu²係数も零となり、従来と同じ `(v,H)=(-3,9),(-4/3,7/3)` に限られる。各枝で残した第1〜4列の六つのminorを取り、t,Pと非零QQ定数だけを除くと、保存された一変数 Bezout 和が1となる。これで特殊枝も除外される。

source replayは元804項からpairを再構成し、旧全域 proof、今回のunused-col1 minor、inverse、特殊枝の全6minorとBezoutを検証する。GB・追加有限体・新しいパラメータ探索は用いない。

## 2. 五二次式の linear syzygy は、ternary pair gcd gate 一つに縮約する

五式を `(Q0,Q1,Q2prime,Q3,D)` とする。Q1,Q2prime,D は最初の四線形座標 `(A,X,W,B)` だけを含む。残る二つを vertex座標 S,U と書く。

**追加仮定 G2:** ternary pair `Q0,Q3 mod(A,X,W)` が互いに素である。これは現在未証明である。rank2だけとは異なり、二つの非比例 conic が共通lineを持つ場合を除く必要がある。

任意のlinear syzygyをpair空間へ制限すると、G2からQ0,Q3のlinear係数は `(A,X,W)` に入る。従ってそれらの係数K0,K3はB,S,Uを含まない。syzygyのS²,SU,U²係数、およびB*S,B*U係数を比較すると、今回証明した2×5 matrixを係数に持つ `c0*K0+c3*K3=0` を五本得る。そのrank2からK0=K3=0。

残りは3つの2-minorのlinear syzygyである。binaryDの既知unitにより(Q1,Q2prime)はheight2（別CI-open noteのproofは全aPJで成立）であり、Hilbert–Burchによって正確に次の二つが生成する。

`W*D-40X*Q1+A*Q2prime=0`,

`B*D+c*C2*Q1-c*C1*Q2prime=0`, `c=6075`。

よってG2が成立する点では、五式のlinear-syzygy空間は正確に二次元である。

## 3. 原四式の quadratic-cofactor rank78 への補題

さらに **G1:** `Q0,Q3 mod(B,W)` が互いに素、と仮定する。G1の失敗は既存の有限Σ coverに含まれ、その上の別計算が必要である。

四式のquadratic-cofactor syzygyをB=W=0へ制限する。G1によりQ0,Q3の二つの係数は定数倍 `(Q3,-Q0)` となる。対応するKoszul syzygyを引いた後、それらは `(B,W)` に入るので `B*L_i+W*K_i` と書ける。W=0でBを除き、さらにA=X=0へ制限すると、G2から L_i は(A,X)に入る。Q0/Q3とQ1/Q2primeとの四つのKoszul syzygyを引いてB*L_iを消せる。

残る二係数はWで割れ、他の二係数は `40X*L+W*K1`, `-A*L+W*K2` と書ける。全式をWで割ると五二次式のlinear syzygyになる。上の二つのscroll syzygyのうち第一はこの復元で零、第二はQ1/Q2prime間のKoszulとなる。従って元のquadratic-cofactor syzygyは六つのKoszulだけである。

四式は独立なので六つのKoszulも独立。六変数の二次cofactor84本に対する行列は `84-6=78` のrankとなる。これは全CIを証明する必要のない、degree4だけの十分条件である。現時点ではG2未証明、有限Σ上も未検証なので、全域quad78は主張しない。

## 4. 具体的な次の有限証明問題

G2の失敗ならternary conic pencil `det(Q0+sQ3)` の四つの係数はすべて零になる。逆は仮定しない（共通kernelを持つsingular pencilもありうる）。actual pairの全係数から計算した四式を、既知t,P,J因子とQQ定数だけで正規化すると、項数は28,56,91,131、全次数は4,6,8,9、t次数は2,3,4,4である。profileを保存した。今後はこの小さな必要locusの有限性とleading-coefficient caseを扱える。未知因子をunitとして除かない。
