# B18 actual D の線形nil二次元と物理的cotangent基底

2026-11-03。全NイデアルのR曲線への制限が正確に(t⁴)であることを独立検証済みの入力とする。
被約横断環のedimは18、全被約coneのedimは37となる。
N⊂D_redなのでT_N⊂T_Dred、両者の次元37が等しいことから **T_Dred=T_N** である。
以前の単なるT_NとT_Rの和からの推測ではなく、DVR fiber productと正確な接触を用いた結論である。

既受理の実150×39 D接平面と実150×37 N接平面を読み、共通39個の物理係数でNを表示する。
全150行の包含を確認した後、その転置のkernelを有理数上で正確に求めると

    u46 + u95,   u63

の二つが基底となる。添字はcanonical bivectorの0-based raw150順序である。
具体的には第一がP^{hw}のhx係数とP^{xw}のx²係数の和、第二がP^{xy}のhz係数である。
P0でいずれも零。scalar chart u106=1でも同じ二つのcotangent classを使える。

一般の還元A→A_redに対し、そのcotangent kernelは

    (Nil(A)+m_A²)/m_A²

である。従って上の二classは、実Dのnilradicalの一次可視部分の正確な基底である。
各classには形式的nilpotent元へのliftが存在する。
**生の二つの線形函数そのものがnilpotentだとは主張しない。**
必要な高次補正、nilradical全体の最小生成系、冪零次数、全Dのdepthは未決である。

再現は既設環境で `python code/ambient_B18_nilcotangent_basis_20261103_v63.py --output <new-dir>`。
源の高次jetやCASは再実行せず、既に独立受理された全接平面とk=4を結び付ける短いQQ線形代数である。
