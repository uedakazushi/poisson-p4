# 全22次不変式の除外と次数下界24：独立監査

2026-11-03。対象 `Q_NEW_FAMILY_DEGREE22_INVARIANT_EXCLUSION_20261103_v63`。
既存36個の QQ source 点の記録だけを使い、選択済み21点で独立な exact QQ 検査を行う。
新しい source 点・beta 点・素数・degree24 probe は一つも生成しない。

## 1. 候補の普遍性と次数

候補24式の定義を旧基底記録から独立に照合した。九つは F4 と既受理の Inv18 trace9基底の積、
十五個は Inv8の最初の五trace4基底と Inv14の三基底の積である。
Inv18の選択は necklace9の列 [6,10,12,18,19,21,23,26,32] と一致する。
Inv14は正確に bracket pairs (1,2),(1,3),(2,3)であり、旧十bracket列では [4,5,7]である。
長さ7のtraceと置き換えていない。

theta,Hの同時共役性、Omega保存、constant Poisson bracketの不変性は既存独立監査を依存とする。
従って trace4は8次、trace9は18次、F4は4次、二つの8次式の Poisson bracketは14次で、
24候補はすべて周囲W上の普遍的な22次Sp4不変式である。
点に合わせた式の発明や、NEW25上だけの関係を周囲の不変式とする操作はない。

## 2. 古い点と全 native 入力を検査した

元36 QQ記録は受理済み degree20 packet と一字義的に同じであり、probeの36記録も一致する。
今回 exact verifier に渡すのは、その中の選択済み21行だけである。
保存されたsource値の pivot、Pfaffian はすべて非零で、各点の所属には既存の全dense-map監査を使う。
このレビューが source DAG を36回再評価したとは記さない。

Omega、その逆、元の10個のrho、sp4基底から別に作ったtrace Gramの逆、word、bracket pair、
全21 betaを再構成した。author native-inputと全tokenが一致し、隠れた点や基底の入れ替えはない。

## 3. 別のexact勾配算法で全値を再計算した

author native codeは実行・importしていない。別のFLINT有理行列verifierで、同じ定義から
V,mu,theta,Hを新たに計算した。微分では dmu を二つの項の直接和から作り、
H=theta+Omega^-1 V Gram^-1 V^Tをそのまま微分した。

各8次traceの全16方向微分は、cyclic contraction K0,K1を使わず、各因子のforward dual matrix product

`(P,dP)(M,dM)=(PM,dP*M+P*dM)`

で得た。この方法はauthor exact codeのcyclic trace contractionとも、probeのadjoint contractionとも異なる。
三つの勾配から Poisson bracket を作り、全378因子値と全504積値を **QQ値そのもの** で比較した。
全て一致した。mod32003で一致したことだけを独立検査の結論にしていない。
新しい微分値も独立receiptに保存する。

選んだ21×21 minorは独立FLINT determinantに加え、別のPython Fraction Gaussian determinantでも再計算し、
完全な有理数が一致した。分子4922bit・分母3186bitで非零である。
32003での3746という値は付加照合に留まり、QQ非零の根拠ではない。

## 4. 全空間の単射性と下界

全22次不変式空間の次元21は、先の独立 two-adjoint reviewの16weight完全DPとWeyl8項から受理済みである。
今回の21個の制限がQQ上独立なので、24候補のspanは全Inv22となり、制限写像そのものが単射である。
したがって NEW25 first-jet hypersurface 上で消える非零22次不変式は存在しない。

この像がSp4安定・既約・斉次な15次元hypersurfaceであることは既存の全source証明を使う。
その生成方程式はsemi-invariantだが、連結半単純Sp4に非自明characterがないのでinvariantである。
既存の次数下界22と今回の22次除外により、次数22は不可能となる。
さらに中心−IはW上で−1として作用するため、invariantの次数は偶数であり23も不可能である。
従って **定義方程式の次数は24以上** である。

24次式の存在、24次候補の基底や評価、別の24次判別式との同定、全不変式環や全fiberの尽くしは未主張。
今回の検算から次の次数を自動的に探索しない。

入口 `code/q_new_family_degree22_invariant_exclusion_review_20261103_v63.py` と独立native
`code/q_new_family_degree22_review_20261103_v63.cpp`。
既設Cloud activation後、新しいoutputへ実行する。60秒・2GB capで、primary資料は変更しない。
