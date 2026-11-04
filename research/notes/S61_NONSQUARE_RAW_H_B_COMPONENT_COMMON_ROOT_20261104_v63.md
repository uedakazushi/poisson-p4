# 成分function fieldでの共通d因子とFULLU境界の確認

小さな還元入力a⁵H,g,kについて、L上のQQ(e)とQ上の2次体QQ(e)[a]/Qで
extended Euclidを行う。Qの既約性は別のexact因数分解packetに依存する。
各除算のwhole quotient/remainder、非零LCとその逆元、最終3本Bezoutと全sourceの
共通因子によるliteral quotientを保存する。単なるGCDラベルには依存しない。

係数はa次数0/1とeのQQ有理関数として記録する。全ての分母polynomialを保持し、
その消失点をgenericな計算で捨てない。各field-polynomialの総supportと保存量を
制限する。これはglobalなFULLU単元逆行列ではない。

さらに共通因子からgcd(U)を反復除去する。全因子が除ければ、generic共通根は
U=0上だけにある。残る分母例外は有限個のe値。Lはaがeの多項式、Qはaモニック
2次であり、実際のa非零点ではHがdモニック4次なので、その例外は有限点だけ。
成分上の空性・例外が実在するか・FULLU全域rank3はここから主張しない。

両成分で境界だけと証明できた場合に初めて、既存の有限残余と合わせB-bad全体の
有限性が従う。片側が未確認又は失敗なら、その条件は未解決のまま保持する。
