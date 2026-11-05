# T1のC単元valuationは正確に四

2026-11-05。fresh ROOT_SOURCE_CODE_REVIEW後、一回だけT1専用Hasse workerを実行した。
未実行のT、他の因子、他tail、商・GCDを追加計算していない。

U1=T1/a³、C=100(a−a0)、a0=(-18e²+9e+2)/100とする。
保存済みT1のliteral7595整数係数と一個のDから、U1(a0+z,e,d)の係数を順に求めた。
**k=0,1,2,3は全係数零、k=4が665項で非零**となった。
従ってprimary計算はv_C(U1)=v_C(T1)=4を証明する。
first layerは逆順減算で零へ戻してから再現し、全665係数をreadbackした。
これはself-consistencyであり、独立Horner/QQ source監査は別gateとして残す。

実行0.118秒、native inner products609591、peak owned33702、入力込みRAW上界914061 bytes。
固定された35188 cells・518 bits・45秒soft/60秒hard・2GB・RAW19.44MBを維持した。
product counterはnative inner loop scopeであり、全scalar/root操作のscopeはROOT addendumに従う。
部分prefixやguard停止は生じていない。

first-nonzero-Hasse.txtはT1CHV001で、k、global denominator D1*100²⁰、
全(e,d,integer coefficient)を保存する。QQ分数を重複展開しない。
completed-zero-orders.txtは四つの完了したZERO層を保存している。

CはFULLUの既知単元である。このvaluationだけでは不良集合やT=0補集合を除去しない。
商T1/(a³C⁴)はまだ計算していない。次は既知因子による膨張を除くための
固定有限商operatorと全remultiplicationのsource予算を別途準備できる。
