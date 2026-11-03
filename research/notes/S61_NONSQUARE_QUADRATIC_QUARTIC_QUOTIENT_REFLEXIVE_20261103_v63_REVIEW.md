# 四次quadratic-product商の反射性：独立数学監査

2026-11-03。rootの `S61_NONSQUARE_QUADRATIC_QUARTIC_QUOTIENT_REFLEXIVE_20261103_v63` を、凍結済みの二つの独立receiptを依存として受理する。新しい算術計算は行わない。

係数環R=QQ[a,e,d,1/(adFBCJ)]は次元3の正則Noetherian integral環である。全22 syzygy独立レビューは、168→210のkernelがrank22のsubbundleで、抽象像Eがlocally free rank146であることを確認した。従って **0→E→R210→Q→0** は完全である。EをR210内のsubbundleとする主張ではなく、特殊fiberでの像rank146を仮定しない。

有限欠損の独立レビューにより、Qは有限閉点集合を除きfree rank64である。有限集合内の各局所環は正則局所環で次元3、EとR210はfreeなのでdepth lemmaからdepth Q≥2を得る。全体のprojective dimensionも局所自由resolutionから高々1である。

Qのtorsion部分Tは、off-finiteでQがfreeであるため有限集合に支持される。Tが非零ならある閉点で有限長部分加群が存在し、その閉点がQのassociated primeになる。これはdepth Q≥2に反する。従ってQはtorsion-freeである。codimension1,2ではfree、codimension3ではdepth≥2なのでS2を満たす。normal Noetherian domain上の有限加群に対するtorsion-freeかつS2の判定から、Qはrank64のreflexive moduleとなる。

長さ1のprojective resolutionによりExt^i_R(Q,R)=0 (i≥2)。Ext^1はoff-finiteではfree moduleのExtなので零、従って有限閉点集合に支持される。これらはsourceの点別rankを推論せずに得られる。

**限界も正しい。** 三次元正則環のreflexive moduleには孤立した非自由点があり得る。例えばQQ[x,y,z]で列(x,y,z)による写像R→R3のcokernelはrank2、depth2のreflexive moduleで、原点では非自由である。この一般代数の例は実Jacobi moduleを同定するものではなく、反射性からfree/flatを導けないことを示す。

従って実四次商の全点free rank64、有限例外のambient rank146、actual C122の全射、全quartic全射は本補題からは従わない。有限Sigma上の新rank78 packetはこの反射性の証明に必要なく、依存に追加しない。

付随するflat JSONはauthor note/receiptと二つの独立dependency receiptのhashを固定し、上記純数学監査の受理範囲を記録する。
