# Cprime の実三次 initial witness と長い指数2：独立正の監査

2026-11-05。受理済み85本の実一次 syzygy と cached V2 の pure-plane 係数から、各 full initial form の pure-plane 射影を独立に掛ける。folded adjoint との85値の一致を確認し、非零は `[49,57,59,77,84]`、最初の49は `99/64` である。

witness49を元の24本の literal V source 列と九変数へ明示的に戻し、`V1*u1=0` の全135成分を再確認する。V0=0かつu1は homogeneous 一次係数なので、真の三次 initial は V2*u1である。V3・未知 constant-kernel 補正・新 source・gauge を使用しない。nonplane 二次項は残余 cubic projection に消されるため、既存 pure-plane cache だけでその商像を正確に求められる。

旧26固定列にこの実 source を追加した27次行列を作り、Schur block inverse を有理数で構成する。両側27次積を直接確認し、rows=旧26の後に25、columns=旧26の後にwitness49という順序で determinant `(-61/18000)*(99/64)=-671/128000` を得る。新しい rank 探索や逆行列 solve は行わない。

よって実際の initial 関係が Cprime/ell³ の全三次空間を張る。Nakayama から `m³(Cprime/ell³)=0`、実 HF333から長さは9である。既知DVR rank2・torsion最小生成数2・Cprime/ell²長さ7は指数集合 `{1,e>=2}` を与え、`length(Cprime/ell³)=7+min(e,3)=9` より **e=2** が従う。

これは Cprime の torsion 指数 `{1,2}` と有限商長さ9の結論である。元の nilideal の nilpotency index・全 nilgenerator・genuine connecting target の実基底や作用行列はまだ構成していない。
