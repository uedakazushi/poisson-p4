# 残二line＋五fixedの全Ore source診断

2026-11-08。producer exact QQ算術PASS、全G閉鎖は未了。
`CLASSIFICATION_PLANE_NONREDUCED_M3_MOVING_SIMPLE_ORE_20261108.json` が全七sourceを構成した。
旧negative factor-rowspanを繰り返さず、全source Euler/Q/strictJac/actualD3lineを保存する。

六源ではweight零及びpositive Euler gaugesで全vertex依存が消えた。
二whole line、DEG1_TRIPLE_EXACT_B、DEG2_CB2_1はQ coefficient space四次元。
二c-dependent壁はQ space九次元で、Q0はw非依存∂w因子、
strictJacobiはcentral F三係数を三角的に解く三式となる。
fullQ/W十八source map及び必要actualD3 line equationsを全て保存した。
Q/W十八のsupersetから元sourceの相関を捨てた結論は主張しない。

唯一DEG2_CA2_0（Croot=a²、A=c²）はpositive vertex-b Euler後も

`b*(p0 z+p2 w) ∂a∧∂c`

がQに残る。この全coefficientは元fullP1kernelの二free座標であり、零へ固定出来ない。
reportのsource_covered=falseを保持し、c-independent Ore構造をこのsourceへ仮定しない。
weight零Q0はzとI=a²w+c²zをCasimirに持つが、higher sourceにその性質が延長するとは未証明。
R=(p0z+p2w)∂a∧∂c自体はQ0とstrict-Schouten可換で、RのJacobiも零。
したがってこの残vertex termを単にJacobiの最高次数だけで零とは言えない。

次の数学は六Ore源の全Q strict-support/actual-normalizer還元と全W module、
及びCA2源の全b-dependent extensionに分ける。
この分岐は原全source上の実際の阻害であり、chosen sliceの取り直しではない。
全generic triple、m二、reduced Rは引き続き別義務。
