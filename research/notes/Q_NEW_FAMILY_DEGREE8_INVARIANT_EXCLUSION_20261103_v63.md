# 新一次像の全8次不変式を除外する厳密証明

2026-11-03。`R4_POSITIVE_WEIGHT6_25`のprimitive beta像Bは、16次元表現内の既約Sp4安定な斉次hypersurfaceである。前の完全QQ Lie-kernel計算で方程式次数≥8を得た。この追加証明は**全8次候補を除外し、方程式次数≥12**へ進める。方程式そのものはまだ構成していない。

## 不変式空間の完全な次元

既に保存したW16の全16 Cartan weightsから、Sym^d(W16*)の文字を整数dynamic programmingで全て数える。各次数の全weight multiplicitiesの和はbinomial(15+d,d)と一致し、符号反転と交換による全Weyl対称性も検査する。C2のpositive rootsを(2,0),(0,2),(1,1),(1,-1)、rho=(2,1)としたWeyl denominator formulaより、trivial representationの重複度は

    sum_{w in W(C2)} det(w) * multiplicity_d(rho-w*rho)

である。全8個のsigned permutationsと各整数項を保存する。次数0–10で非零はdegree0:1、degree4:1、degree8:5だけである。これは標本によるrank上界ではなく、有限な文字の完全計数とWeyl character theoremである。degree2/4/6は既存の完全QQ Lie-kernel結果0/1/0とも一致する。

## 明示的な8次不変式と一致する下界

不変交代形式omega、pullback anti-representation rho(A)、trace Gramからmoment map muを作る。F=Tr(mu²)/2、theta=rho(mu)、H=omega^-1 Hessian(F)とする。全十Lie作用がomegaを保つこと、全100 commutatorが正しいpullbackの符号を持つことをQQで検査した。momentのequivarianceとFのinvarianceからthetaとHはいずれも同じconjugation作用で変換する。従ってtheta/Hの四文字cyclic trace六個とTr(mu²)^2は、すべて斉次8次不変式である。Hは係数展開を省略しても、有限公式

    H = theta + omega^-1 * V * Gram^-1 * V^T,
    V_A = omega*rho(A)*b

で完全に指定される。

実dense mapの五つの正則QQ source点で、七個の不変式の評価5×7行列がrank5となる。選んだ5列の非零QQ determinant、全beta24、value pivot、Pfaffian、全評価を保存した。よって五個はBへの制限後も線形独立である。全不変式次元が5なので、これらは全8次不変式の基底となり、そのBへの制限は単射である。従ってB上で零となる非零8次不変式は存在しない。

特にmoment determinantとquartic squaredの比も最初の二点で異なり、det(mu)-c*I4²という任意の定数cの候補は除外される。だが結論はこの二次元pencilだけでなく全五次元8次空間に及ぶ。

## 方程式次数と範囲

Bの主prime idealの生成元はSp4のcharacterが自明なため不変であり、scalar安定性から斉次である。中心−Iによりodd degreeはない。degree2/6/10は不変式がなく、degree4は既に非零評価で除外し、degree8は上の単射で除外した。よって最小方程式次数は少なくとも12である。

この結果はBの方程式、全高次数不変式、一般S61 fiberの全点列挙を与えない。五点評価は不変式の独立性と候補除外だけに使い、恒等零やgeneric rank上界を標本から推測していない。

保存先: `code/q_new_family_degree8_invariant_exclusion_20261103_v63.py` と `reports/Q_NEW_FAMILY_DEGREE8_INVARIANT_EXCLUSION_20261103_v63/`。前のrecognition packetは変更しない。
