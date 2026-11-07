# cubic scroll・integral零曲線の排除：独立確認

2026-11-07。rootの幾何提案を `log_manuscript` が別に再導出した後、rootが完成証明を再確認した。
純粋な幾何命題は **PASS**。同note末尾の21源/11変数計算は未閉のproducer結果であり、
このreviewはその数値上界を改善したとの主張ではない。

確認した全矢印は、rank二conormal可換、全nonplane位数二、Hamiltonian-parallel G/sigma、
全primitive rootとvertical零のfactorization、rank一のC0係数parity、rank二のV0/R両枝である。
V0ではE上の一般ample complete-intersection curvesを使い、全codimension二bad setを避ける。
数値Gysin関係だけでA1−A2のH及びH+fとのintersectionが零となり、Pic(F1)の二座標を決定する。
Rの全点smoothnessやΓの全点SNCを使っていない。

V=Rではdet rootのglobal nowhere-zero sectionから正規束の全域splitが必要。
二つのnormal exact sequenceをparameter上で確認した。
rulingのO→O(1)^3は二つの独立linear formを使い、quotientはO(1)+O(2)。
directrixのO(-1)→O(1)^3は(s²,st,t²)で、二つのlinear syzygyからquotient O(2)^2。
したがってsplit候補は(C0+3f)+(2C0+4f)だけだが、積は八、実c2は自己交差九。
この矛盾で全V=Rを排除した。

normal quadraticのrank一/二と全vertical場合は尽くされ、integral reduced Rの全枝が空。
Rの可約非被約、cubic cone、全rank二scroll/Gの尽くしは残る。
二回の未閉linear/leading boundを成功理由に使わず、余分なGB又はjetも要求しない。
