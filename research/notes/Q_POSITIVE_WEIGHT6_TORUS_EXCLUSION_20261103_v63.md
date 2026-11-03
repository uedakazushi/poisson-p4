# 正の重み六軸はT02・T03にも入らない

2026年11月3日。全eta非零の正軸のprojective PGL5軌道次元は既証の22である。
一方、既存の全MPS source表示から、T01の全点の軌道次元は20以下、
T02・T03の全点の軌道次元は21以下である。従って正軸は三成分の全閉包を避ける。
T01については全274関係による別の排除と一致し、T02・T03の排除が新しい。

T02・T03のdenseなfull sourceは対数型tensorに一つのedge smoothingを加えたものである。
edgeは(0,1)、余りの指数はそれぞれ(0,0,2)、(0,1,1)。対角PGL torusの次元は4で、
edgeのcharacterはそれぞれ(1,1,0,0,-2)、(1,1,0,-1,-1)である。
対数項は全対角torusで不変であり、smoothing項もそのcharacterのkernel上で不変である。
したがって両sourceの全admissible tensorは少なくとも三次元torusを安定化群に持つ。
T01ではedgeがなく、四次元torusが残る。

これは一つの標本のstabilizer計算ではない。一般periodと全edge変数を保つ
full source上の恒等的な対角actionであり、そのPGL飽和像は対応成分に稠密である。
無限小actionのrank上界は全小行列式の消失という閉条件なので、その成分の全境界にも延びる。
正軸のrank23のGL+scalar小行列式はこの閉条件に反し、T02/T03の非所属を証明する。

T04以降では残るtorusの次元が二以下なので、この同じrank上界による排除は得られない。
他成分の帰属、新成分、全局所枝の尽くしを主張しない。
