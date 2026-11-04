# Qの分母例外を保護した15normによる有限cover縮約

rootのnormal formはgeneric fieldでの式なので、分母例外上で全root-normを必要条件
とはみなせない。まずexc=gcd(SF39,productE_Q*D)を全て残す。残りの互いに素な
部分だけを15個の非零minor-root-normとGCDし、excとの積を新しい安全なcoverにする。
盲目的にSF39全体と15normのGCDを取る操作は行わない。

全GCDについてBezoutと両方のliteral quotientを保存し、独立整数演算で確認する。
新coverは旧SF39を割り、全actualB-bad候補を保持する。これで有限代数の次元上限を
小さくできるが、fiberでのunit検査やB-bad空性はまだ証明しない。
新たな点・素数・GB・nested resultantは使わない。
