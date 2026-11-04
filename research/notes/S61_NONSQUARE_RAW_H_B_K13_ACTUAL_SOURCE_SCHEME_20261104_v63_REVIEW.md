# Actual I3(B) の Ka13 同型を与える両写像の独立監査

対象は凍結済み ACTUAL_SOURCE_SCHEME。作者の Ring・worker を呼ばず、別実装の
Ka13 固定 monomial-table 演算と通常の整数多項式恒等式で確認する。新しい逆元・GCD・
solve・GB・点・素数は探索しない。

全312 monomial nodes の literal inputs と積の整数 g13 quotient を検査し、monomial
map の全 exponent labels を追跡する。元の全22 source generators（2963項）を、還元済み
次数12以下の係数の線形和として別計算し、厳密に零と確認する。FULLU の49元項の評価も
同じ元 source から確認する。

旧 NAMED_SATURATION の七つの整数 Bézout・content・inverse products を再検証する。
四つの parameter unit は元 a,F,C,b1 の式に一致することも確かめる。r=0,rB,rJ について
H4(d)=(d-r)Q_r(d)+H4(r) を全係数で検証し、H4(s)=0から保存された(s-r)逆元を再構成する。
FULLU=a*d*F*B*C*J の全六因子、先頭単元、最終逆元との積1も確認する。G(a,t)の49 source
項と既存整数 Bézout は別に検査する。

scheme の upper map は、G(a,t) の新しい逆元から推測しない。既に独立受理された全G0
除外が G を R/I で単元にする。受理済み norm source identities から H,h1,k1,m1∈I、
exact Ka13 source bridge から g13=0,e=t、whole-Ka13 d-Euclid から d=s が従う。ここで
R=QQ[a,e,d,1/FULLU]、I=actual I3(B)で、22 normalized generatorsとの ideal 等号は
既存35minor source reviewを保持する。今回の明示 U逆元は局所化の生成元1/Uも a の
多項式で表す。従って alpha:Ka13→R/I は well-defined な全射である。

逆向きの beta:R/I→Ka13 は、今回の全22零性と FULLU逆元だけで well-defined な単位的
評価写像として構成される。beta∘alpha は生成元aを固定するので恒等。従って alpha は
同型であり、actual B Fitting bad algebra は正確に QQ[a]/g13、QQ長さ13で非空である。
体性、既約性、reducedness、étale性、幾何点数は仮定も結論もしない。

対象は actual7×3 B の Fitting bad scheme に限る。coker B の module length、full
cubic/Poisson rank loss、quartic239、generic exhaustion へは移さない。入力 hash・実行結果は
report と independent-proof に保存する。
