# rank零plane・normalnil非可換affine型は全重複度でG非優越

2026-11-07、主義務G。producer幾何の追加橋。旧m二certificateは変更しない。

**命題。** actual退化因子が被約既約quinticで、rank零plane上のtraceless first symbolが
非零global vector field Vとconstant nilpotent normal sl2 matrixの積であるとする。
trace field TとVのweight零Jacobiが`[T,V]=−bV`、b≠0の場合、
全higher Poisson sourceのambient projective像はdimension二十二以下。
planeに沿うDの重複度二を仮定せず、三・四も保持する。

`REDUCED_DEGENERACY_ZERO_SURFACE_ABELIAN_CONORMAL.md`とその別作者reviewにより、
被約Dでは全重複度でconormal可換、従ってfull35/30/9 sourceを使える。
`PLANE_RANK_ZERO_GLOBAL_VECTOR_NILPOTENT_FIRST_JETS.md` §§1–2の全Jacobiと
actual shear/scalar/centralizer normalizationはm二を使わず有効。
Kのaをnormal shearで消しb=1へ直すと、VはJordan3又はJordan21。

Jordan3のwhole sourceはV=J,T=diag(1,0,−1),K=zwへ直る。
そのleading determinantは常に非零なので全てm二であり、既存full-weight certificateで二十二以下。

Jordan21の全sourceは

`V=N=E01`, `T=Dρ+xN+yB+zC`, `K=zw`,

`Dρ=diag(ρ+1,ρ,−1−2ρ)`, `B=E02`, `C=E21`

で、Eijは係数行列、対応fieldはxj∂xiである。
全leading determinantは

`det(Euler_base,T,V)=x1²*(z*x1−(1+3ρ)*x2)`。

従ってm三以上となる全sourceは正確に`ρ=−1/3,z=0`。
この場合B方向は非共鳴なのでactual exp(sB)でyを消し、exp(qN)でxを消せる。
逆もfinite unipotent change。whole sourceは

`T=D(−1/3)`, `V=N`, `K=zw`

という一つのnormal formへ直る。
選んだsliceだけでなく、全rank一nil V orbit四、normalnil line一、Kのa,b二、
残るT解二（x,y）を復元したfixed-plane affine symbol sourceはdimension九。
これは全source十一からρ条件一、resonant row零条件一を引いたものでもある。

immutable report `CLASSIFICATION_PLANE_RANK_ZERO_GLOBAL_LOWPARAMETER_WEIGHTS_20261107.json`の
`NIL_B1_DIAGONAL`、rootρ=−1/3は、その当時m二外として捨てずにfull operator算術を保存している。
全weight一D1 kernel七、全weight二D2 kernel零である。
新modular入口により、これを今回の全m三以上sourceへ適用出来る。
各nonempty higher fiberはdimension七以下なので

`dim image≤symbol9+higher7+plane6−scale1=21`。

他の全m二branchesは既存whole低parameter証明で二十二以下。
ゆえに本affine型は全重複度で二十二以下となる。

本命題はnormalnil b零のcommuting/pointwise-dependent source、normalss highm、
scalar traceless symbol零、又はk一のhigher-rank leading-Pf零を含めない。
全plane又は全Gの完了ではなく、旧表で除外した一つのwhole補集合を正当に復活させた結果である。
