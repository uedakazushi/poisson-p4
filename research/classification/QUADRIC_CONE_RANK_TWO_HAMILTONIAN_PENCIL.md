# rank 二 cone の正しい Hamiltonian pencil と全 Jacobi 条件

2026-11-07。主義務 G、producer `plane_collision_review`。
`QUADRIC_CONE_RANK_TWO_TANGENT_SOURCE.md` の全69raw sourceを出発点とする。
初回の有限実行はsemisimple Cの全modh affine source（Euler quotient後33次元）まで成功した後、
vertex変数に依らない旧rankzero Ore形を仮定するassertで停止した。
成功とはせず、失敗時のproducer bytes・全33方向を
`CLASSIFICATION_QUADRIC_CONE_RANK_TWO_SOURCE_RUNTIME_20261107.json` に不変保存した。
初回script・jobも変更しない。

## 1. 診断：旧 Ore identity は全 tangent source の identity ではない

rank二restrictionの全F quadraticから \(F_v=\ell(u)\) が残る。
保存済み全kernelでは、C=∂0∧∂2のとき

\[
 B(u,v)=B_0(u)+vD,\qquad
 D_{01}=-2\ell_2,\quad D_{02}=0,\quad D_{12}=2\ell_0.
\]

例えばkernel内のF=u0v方向はD12=2を持つ。
これは全modh必要sourceに対する反例で、full Poisson solutionの存在例とは主張しない。
full Jacobiが後でD0を強制する可能性は残すが、入口でD0としてsourceを切ってはならない。
同様に \(\operatorname{Jac}(B_0)=0\) を単独で課すことは不正である。

## 2. 正しい全 affine 形

chart h1の全sourceについて検証すべき形は

\[
 P=B_0(u)+vD+(Y(u)+vT(u))\wedge\partial_v, \tag{1}
\]

Dは**constant bivector**、\(\deg B_0\le2\)、\(\deg Y\le3\)、\(\deg T\le1\)。
この全33source identityを、semisimple/nilpotent/zeroの全C型で確認する。
旧rankzero OreはD0の特別枝に当たる。
Fvv零、全Fv三係数、Fu-only五係数も完全kernelから復元する。

## 3. 正しい五つの全 Jacobi block

constant DゆえJacD0。
(1)の全affine Jacobiをvの係数で展開すると、正確に

\[
\begin{aligned}
 \operatorname{Jac}(B_0)+Y\wedge D&=0,\\
 [B_0,D]+T\wedge D&=0,\\
 \mathcal L_YB_0-Y\wedge T&=0,\\
 \mathcal L_TB_0+\mathcal L_YD&=0,\\
 \mathcal L_TD&=0.
\end{aligned} \tag{2}
\]

最初の二つはbase三方向のJacobiのv0,v1、残りはmixed三方向のv0,v1,v2。
全sourceをpolarizeしてこれらの全係数とoriginal350 projective Jacobiを比較する。
quad(a,b,c,h)のh1 restrictionの負号がaffine Jacobiであり、global Euler relationとdense chartにより同値。
二項目を省いて旧Oreの三blockだけにすることはない。
full126 actual-Pfaffianも同じ全source上で保持する。

## 4. 数学的 route change と未解決の範囲

初回失敗のminorやjetを増やす操作ではなく、全sourceのHamiltonian pencil構造への修正である。
新しい一回の有限source実行は同じ全59unknown、全conormal十一式、全350modh係数、Euler五を保つ。
affine kernelを狭める誤った \(\operatorname{Jac}(B_0)=0\) は課さず、(2)と全original equationsを保存する。
ここでの成功は全source入口の完成であり、G≤24を得たとは主張しない。
次は(2)のwhole support/linear-fiber bound又は明示的既知family認識が必要。
GB、未知parameter点列、normal jetの増加は行わない。
