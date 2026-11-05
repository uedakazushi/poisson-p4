# Cprime の既存 initial による三次上界 gate

完全な二次作用素の階数 1 により、Cprime の真の Hilbert 関数は次数 2 まで (3,3,3) である。本 packet は三次の値や階数を計算せず、次の小さい正の上界試験の入力を固定する。

二次 NF4 の階数 1 の像を消去して NF3 を得る。実際の二次 kernel 列 18 の pure-plane 部分 P は、同列の完全な二次 initial と NF4 が一致する。NF4 は pure-plane 以外を読まず、その核は既知の一次式の一次倍の空間に等しいため、P 自身も真の二次 initial relation である。NF3 の核はその既知空間と P の和である。

三次自由空間の各単項式を、ソートした最初の二変数を二次 divisor として NF3 に送る。標準三基底の九変数倍を target とし、I-qE の全 27 consistency 列を残す。44 本の一次式の二次倍 1980 列、P の一次倍 9 列、および三 target の ell³ を加えた計 2019 列の cokernel は、Cprime/ell³ の既知三次 initial quotient に等しい。

三次 quotient がゼロと証明されれば、Nakayama により m³(Cprime/ell³)=0、長さの上界は 3+3+3=9 となる。一方、既知の DVR rank 2 と torsion 指数 {1,e≥2} から長さは 7+min(3,e)≥9 である。従ってこの正の試験が通る場合に限り e=2 が従う。既知三次 quotient が非零でも e>2 とは結論しない。

将来の worker は二 pass の streaming とし、全 2019 列を保存しない。独立な元列は高々 27 本であり、固定小行列の両側逆・全列 factorization hash・元 source/projection hash を保存する。全入力と出力は 32 MiB、owned QQ は 24000、係数 2048 bit、60 秒、2 GiB を上限とする。完全な二次値の独立 peer と、この別 SOURCE gate の peer および親の GO 前には三次値 worker を実行しない。新 source・gauge・jet は不要である。
