# moving-double Ore hand leading係数の訂正

2026-11-07。旧 `PLANE_NONREDUCED_M3_MOVING_DOUBLE_ORE_REDUCTION.md` のbytesは保持する。
§2のD3 identityのsignを訂正する。元のactual Pfaffian conventionで

`D3=+(3/2)u A(A−1) z(bz+aw)^2`

が正しい。旧手noteのminus signは誤り。
R3Lのexact coefficientは `D3=6u z(bz+aw)^2` である。
旧note§3はnonzero constantとしか述べていないが、最初のcodeの予測−3uも誤りだった。

最初のbounded jobはこの不一致で停止し、runtime failure receiptを保存した。
元P0、Q全source、W、全Pfaffian式は変更せず、予測identityだけを正しい値へ直す。
D3非零条件u*A*(A−1)非零、二excluded parameters、central coboundary、actual normalizer、
全三canonical moduleリストは変わらない。
失敗runを成功証明とは呼ばず、後続成功reportの実際の計算と区別する。
