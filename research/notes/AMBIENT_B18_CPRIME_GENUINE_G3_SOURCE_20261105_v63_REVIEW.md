# cached 実九次元作用と真正 G3：独立 SOURCE 監査

2026-11-05。実作用・kernel値・G3作用は計算しない。same Gamma frame の固定V1二十四列・全九一次行・逆・一次NFを独立に再構成し、旧二次NFと基底を固定する。

M=Cprime/m³では実関係 `V_j=L_j+Q_j` とその一次倍を使う。既知131次元の一次倍と実二次initialの全kernelがkerpi3だから、Qのnonplane部分を除くことは実関係を引く操作である。未知を零とはしない。一次入力 e から実pivot関係を `I24*e_R` 倍だけ引けば、NFは `(0,q1e,−pi3Q_piv*I24*e_R)`。二次入力は `(0,0,pi3f)`。全四十四関係と可換性・cyclic九基底を値workerで確認すれば、実HF333の長さ9との一致により同型になる。

基底順は三head、F=(gamma0ell,gamma2r,gamma2ell)、Q=(gamma0ell²,gamma2rell,gamma2ell²)。ell作用EはF_iをQ_iへ、Qを零へ送り、headの一次blockは(F0,0,F2)である。残る二次blockをSと書けば、固定5次minorはdiagonal1の下三角、E²の二本はQ0,Q2のI2となる。`K=H1−sum S_i1 F_i` を用いてkerEは(K,Q0,Q1,Q2)、imE²は(Q0,Q2)、genuine basisは(K,Q1)、kernel上の射影はambient座標(1,7)である。新しい rank/inverse 探索を要しない。

既証指数(1,2)とM=Cprime/ell³により、このquotientは真正Cprime[ell]である。可換性と三層supportから、各変数はKをalpha_i Q1へ、Q1を零へ送る。全alphaはまだ未計算であり、genuine connecting値も計算しない。

入力の各typed積のdenominatorと絶対値だけを検査し、作用を実際に足し合わせずcommon denominator11547360000000000と181bit上界を確認する。18000 numeric envelope・hard24000・60秒・2GB・32MiBを受理し、未実行workerの独立`--help`だけを実施する。元のnilidealのgenerator・cube・H0長は対象外である。
