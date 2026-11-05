# source に結び付いた七次元 D 商と真正四次元 target：独立監査

主 SOURCE `a55309b4...`、worker `fe20a0bf...`、gate `0bea5314...` に対する
独立報告 `bd84ccc8...` を固定した。四つの response はまだ評価していない。

C4 の実23基底を、A、B、rA、rB の四本の ell 鎖16個と
(H,rH,ellH,r²A,r²B,r²ellA,r²ellB) に独立に組み直した。
次数零〜三の既証 literal NF と、四次の非零 cut minor の補行1、4を使うと、
各 graded block が実基底となり、m⁵C4=0 から23個は実 filtered 基底である。

四本の E³像は、次数三の ell³A/B と次数四の r ell³A/B で独立。
ell³Tor_D C=0 なので、C/Tor_D C から ell³C への ell³ 写像は同型である。
従って四 source は自由商の mod ell 基底となり、D-Nakayama からその全 Dspan F は
primitive D-free 四となる。選んだこの F について C=F⊕Tor_D C と分解する。
ただし D 加群としてだけである。r(rA)=r²A は商の非零座標を持つため、F は R0 安定でない。

q:C4→C4/F16 の七次元商は、実 Tor_D C と D-linear に同定できる。
全 m⁴C4 は F16 に入り、z を含む三次と H の平面三次も m⁴へ上がるので、
q(ell*rH)、q(ell²H) は二つの三次 source -P1³、-P2³ の射影だけで求まる。
curve の r³還元と q の F 消去により、必要な四演算子は旧 C2 の24葉・60項だけとなる。
lambda、mu、delta4 又は新しい jet を零と仮定したのではなく、その座標が商で消える。

四係数を独立形式変数とする E7 の rank三 unit minor、全28 kernel 積、
自由 I4、kernel の完全性と E³=0 を確認した。
q は ker E23/im E³ から ker E7 への同型を与える。
kernel は F16[ell]=ell³F16 だけなので、有限 kernel をそのまま真正 target とする誤りはない。

値 worker の old C2 leaf照合、全 add/mul 前の guard、E*K と E²/E³ の非零 path が
定数倍一次・最大二項だけであることを確認した。全高は43bit、演算前87bit、
静的25000・hard30000・全体150000、RAW/時間/memory cap 内である。
本 PASS は実行準備の前提を閉じるだけで、親の別途一回実行指示が必要な段階を保持した。

商への素朴な r/z 作用、R0 transport、実形式 torsion lift、connecting 値は未構成。
再現は `code/ambient_B18_C_genuine_T7_source_review_20261105_v63.py`。
