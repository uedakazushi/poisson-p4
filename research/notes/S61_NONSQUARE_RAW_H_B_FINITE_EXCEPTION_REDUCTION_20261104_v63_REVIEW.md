# 有限例外の縮約とL全成分の単位証明：独立監査

2026-11-04。対象は `EXCEPTION_SQUAREFREE`、`Q_EXCEPTION_NORM_FILTER`、`L_FINITE_UNIT_FIXED` の凍結三packet。既受理の全B-rankbad有限性・元minor/source gateへhashを結び、著者solver、GCD、点、素数、GBは再実行しない。

## Squarefree supportと保護したQ例外

L/Qの元例外積の次数は65/361。全source factorのQQ単位によるmonic化、全積、微分、保存GCDの両割り切れとBézoutを独立QQ[e]演算で検査する。SとS'のBézout=1、およびPからSの因子だけを除いて最後に1になる全29/59段の掛け戻しを確認する。従って根集合を変えずsquarefree9/39へ縮む。新しいfactorizationやGCDの出力を信じるだけの検査ではない。

Q側ではE=product(19分母)*Dをsourceから再構成する。保存excについてexc*generic=SF39、exc|E、Bézout(exc;SF39,E)を確認し、excは正確に37次である。これら37根は全て残す。互いに素なgeneric二次部分だけを元十五normへ順次掛け合わせたGCD証明の全quotient/Bézoutで確認し、最後が1になることを検算する。従って二つのgeneric根だけを除ける。分母例外上でgeneric root式やnorm条件を無条件に使っていない。

通常QQ[e]の固定恒等式は独立native polynomial interpreterで検査する。逆元・GCD・Euclidは実行せず、係数の積和とexact equalityだけを使う。旧h1/k1のe最高係数が非零QQ定数であることと、そのmonic正規化全係数も確認する。

## Lの有限QQ代数と全22source

S_L9から除くのはa_L(e)=−9e²/50−3e/50+3/25の零成分だけである。保存clopen因子の両掛け戻しを確認し、残るS7とa_Lのglobal Bézout=1、有限代数内のa_L逆元を検査する。a=0はFULLU外なので、この除去でactual候補は失われない。元S9はsquarefreeであり、残りS7上のa逆元も確認するため、未確認のa=0部分を単位扱いしたものではない。

別のQQ有限代数実装で、元a^5Hをa=a_Lへ直接代入し、a_L^-5を掛けてdモニック四次Hを構成する。全五係数を保存値へ比較する。E=QQ[e]/S7、A=E[d]/Hは28次元で、基底はe^i d^j、index=7j+iである。積はeとdのmonic reductionから定義し、著者の行列solverやalgebra helperを使用しない。

元の全22minorとFULLU product Uを直接この代数へ写し、全source imageを照合する。全22×28=616 ideal列を元image×基底から再構成し、17222個の非零係数も確認する。Uを28回掛けた値と、元22sourceに対する保存cofactorのdot productを別に計算し、

`sum_i C_i P_i = U^28`

を確認する。これは点での単元推測でも、matrix rankだけの主張でもない。

Uを反転すると右辺が単位になるため、この有限L coverのFULLU部分ではminor idealは単位idealである。既存のL generic部分がB=0境界にあることと、全S9例外をこの証明で処理したことを合わせ、**L成分上のactual FULLU B-bad点は全て除去される**。初回のcompile-only失敗は元ファイルのままhash-boundで保持し、成功に読み替えない。

## 残る有限目標

Qの保護した37根は未処理である。Qはaモニック2次、a≠0ではHはdモニック4次なので、Q側の幾何学的候補は高々296点となる。旧h/G,k/G由来の高々1600点を保持し、全FULLU B-bad点数の安全な上界は**1896**へ下がる。

これはscheme lengthの上界ではなく、Qまたは旧残余の空性も示さない。全域B rank3、p=0のraw H7、full cubic57、quartic239、Poisson解尽くしへの昇格はない。全input/evidence hash、固定恒等式program、実測値は同名REVIEWに保存する。実行バイナリはignored workだけに置く。
