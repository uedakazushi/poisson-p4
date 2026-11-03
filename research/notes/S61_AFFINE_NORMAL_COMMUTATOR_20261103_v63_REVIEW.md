# Normal-only commutatorと残余表示：独立監査

2026-11-03。全129×126 cross入力から、別実装のSymPy DomainMatrix GF(32003)逆行列で三つの42×42線形係数行列を再構築した。exact QQ native出力はstream hashを照合し、必要な全5292線形係数を別に有限体化した。三つのcommutator行列と独立に計算した行列式は21866,23983,19574となり全て非零である。QQ p-integralityと元134式へのbindingは有限5生成packetの全証明に依存する。

実商環で `z_i*n=R_i(n)n`。係数はnだけの可換多項式なので、二つのzを交換すれば `[R_i,R_j]n=0` が元ideal内の厳密な三次関係となる。modpで非零の定数commutatorはQQ上も非零で、originにおける42方向Jacobianは可逆である。

2796生成元のnormal-only Jは全origin外chartで実残余を与える。42三次式を加えたJ0はoriginでも局所ideal mとなる。残余は有限5生成とglobal splitterから閉かつoriginを避けるので、全極大idealでの局所化を比較して `J0=m∩I_B` を得る。全unknown scheme構造を保存した等号である。

adjugateにより `det[R0,R1]*n=0`。定数detで割ったD84はknown五点上1、nがunit idealを生成する残余上0なので、既存のglobal idempotent eと等しい。従って `I_B=J0+(D84)` がunsaturatedなnormal-only残余表示となる。

一般beta上の有限性、残余unit、空性、全五点の尽くしは結論しない。行列式は次数84以下の有限recipeであり、巨大QQ展開や次数84の等号を計算した意味ではない。元のdegree3 Macaulay障害も再試行していない。
