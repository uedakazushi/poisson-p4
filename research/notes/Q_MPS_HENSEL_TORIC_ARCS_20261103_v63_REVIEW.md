# 十六の正の Hensel MPS 証人の独立監査

2026-11-03。Cloud 復旧後、停止前の fresh report と全証明書の実読出し・一致を確認した。
元 period basis の完全 annihilator と五行列の直接置換から全五整数方程式を独立に再構築し、
十六の Fp 点、全非 rho Jacobian と pure integer Gaussian elimination による determinant を検証した。
全固有値差・rho=0,±1の回避、period 分母の p 単元性も確認した。
独立した整数多項式評価で、固定 rho の全 p^8 Newton 値と各更新を掛け直して一致した。

正の全五方程式の解と可逆微分により多変数 Hensel が真の Zp 解を与え、同じ微分は
rho への source 射が étale であることを示す。rho=rho0+u の形式解も得られるので、
QQ[rho] の非零多項式が admissible source 上で恒等零になることはない。
従って QQ 上で優越的であり、忠実平坦な C への体拡大後にも rho を支配する幾何成分がある。
単元条件は同じ非空開で保たれる。実 period subfamily と actual canonical arc の Qa 極限、
さらに閉性から、全有限複素 a の指定十六成分への所属が従う。

これは有限体の空性や長さを標数零へ外挿する議論ではない。
八桁の近似自体を全次数の証明とせず、all-order existence は Hensel と étale 性から得る。
有理数解を発見したとも、全通過成分・局所枝を尽くしたとも主張しない。
監査結果: `reports/Q_MPS_HENSEL_TORIC_ARCS_20261103_v63_REVIEW.json`。
