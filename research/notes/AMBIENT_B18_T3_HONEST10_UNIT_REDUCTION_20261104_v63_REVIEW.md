# Honest10単位Schur消去の独立監査

2026-11-04。既受理のhonest10 presentationを前提に、全scalar originをtyped DAG semanticsで
再評価し、元10×2070 relation matrixのrank1、top追加後rank1、top12原点零を確認した。
旧regular quotientは受理済みの原点値を使い、792個の新しいcoordinate quotientの原点値は
honest10独立reviewでactual source ADにより既に検証されている。DAG keyの全域一意性は仮定しない。

使用relationは1389列、free generator5の係数が10/3。非零原点を持つため真の局所unitである。
全2069残余関係とtop12、9×10射影の各entryについて

    v_k -> v_k−c_k*u^(-1)*v_5

を作者とは異なる積のmultiset正規化で照合した。これは有限jetのrank消去ではなく、
unit relationでgeneratorを解く正確なmodule同型。全36 Koszul kernelを含む元のpresentationから
移したので、二行ambientへの埋込みでTorを捨てる操作ではない。

全残余relationとtopの原点係数は零。literal zero以外の関係を捨てておらず、今回は零列削除も無い。

    D3=coker(R0^2069 -> R0^9),
    D4=coker(R0^2081 -> R0^9),
    mu(D3)=mu(D4)=9.

このminimal generator数はcokernelのresidue-field次元から正確に従う。
T3はD3内の12列の像であり、そのrelation kernel、非零性、actual membership、depthは未判定。
新しいjet、point、GB、nonunit inverse、作者producer再実行は行わない。
