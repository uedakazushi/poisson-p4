# 新一次像の全12次・14次不変式を除外する有限証明

2026-11-03。前の8次packetを変更せずに続ける。新成分のprimitive beta像Bの斉次方程式は**次数≥16**である。方程式そのものは未構成であり、一般betaのS61全点列挙を主張しない。

## 完全な不変式次元

全16 weightsからSym^d(W16*)の全weight multiplicitiesを整数計算し、C2 Weyl character formulaをdegree14まで適用する。各次数の全weight数はbinomial(15+d,d)に一致する。非零のtrivial multiplicityは0:1、4:1、8:5、12:9、14:3であり、他は零である。全8 Weyl項を各次数で保存する。

## 12次の全空間

前のpacketで厳密にSp4-equivariantと示した二つの2次行列theta、Hを用いる。それらの六文字cyclic trace14個は12次不変式である。実dense mapの正則な九つのQQ source点で評価すると9×14行列のrankは9であり、非零の9×9 determinantを保存する。全12次不変式次元が9だから、選んだ九個は全空間の基底であり、Bへの制限が単射である。したがってB上で零となる非零12次不変式は存在しない。

## 14次の全空間

8次packetの五個の基底をf_iとする。omegaが不変な定数交代形式なので

    {f_i,f_j}=grad(f_i)^T*omega^-1*grad(f_j)

は14次Sp4不変式である。十個のbracketを同じ実QQ点で評価し、最初の三点の3×10行列がrank3となる。全14次空間の次元も3だから、選んだ三個のbracketは全空間の基底であり、Bへの制限も単射である。

微分は有限なquadratic matrix公式からcyclic trace product ruleを使って正確に計算する。degree8関数の一方向微分を−4,…,4の九点のLagrange補間でも独立に計算し、五個全てが解析的gradientと一致することを検査した。これは有限差分近似ではなく、次数≤8の多項式に対する正確なQQ恒等式である。

## 方程式次数

既約Sp4安定hypersurfaceBの生成方程式はSp4のcharacterが自明なため不変、scalar安定性から斉次である。中心−Iによりodd次数はない。前の完全kernel/基底証明で2,4,6,8,10を除外し、今回は全12/14次を除外した。したがって次数は少なくとも16である。

評価点は独立性と候補除外にのみ用いる。標本rankから不変式空間の上界やpullbackの恒等零を推測していない。上界は完全な文字計算と標準Weyl character theorem、下界と制限の単射性は実点の正確なminorで得る。これ以上の次数探索は、このpacketの独立review後に行う。

保存: `code/q_new_family_degree12_14_invariant_exclusion_20261103_v63.py` と `reports/Q_NEW_FAMILY_DEGREE12_14_INVARIANT_EXCLUSION_20261103_v63/`。前のdegree8/低次数packetは不変。
