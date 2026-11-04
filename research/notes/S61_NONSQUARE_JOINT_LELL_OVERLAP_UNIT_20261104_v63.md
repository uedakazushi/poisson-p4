# 四monic二次の overlap filtration と残る二次二本

full183行列を展開せず、既存joint frameの四monic二次だけで正確な三次S行と
四次18 overlapを reduceする。全 Z/W 補正を保持する固定 triangular reduction
であり、GBや単項式idealの同定は使わない。

18x141だけについて既知単元 a,d,F,B,C,J の係数を探し、単元pivotのみで行消去する。
未知因子を割らない。全行変換と全2538係数の等式を検証する。
残る行が全て零なら、そのpivot数だけ全点で追加の image方向を得る。
これは source-kernel方向の追加とは異なるので、既知38・残り42を変えない。

次にq0,q5の29倍行をそれぞれ同じ98 pivotと追加pivotで reduceし、最大58行の
正確な小blockを保存する。ここでは追加のGaussianをせず、個別の既知単元係数数
だけを数える。全112 cubic倍行と全33 genuine quarticは元の正確なlazy oracleに
保持され、この計算で落とす式はない。

一回60秒・2GB・32MiB。全183x141行列、新しい点・標数・GB・generic rankを
計算しない。新たなtarget縮約を使うには、得られたunit行演算をfull sourceへ
有限circuitとして輸送する必要があり、full239は依然未解決である。
