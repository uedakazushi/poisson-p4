# W1 tailのLaurent操作：w多項式algebraを使い、w=0完備環を代入しない

2026-11-04。凍結済みW1 coordinate preflightおよびtail bindingノートの
`D_perp[w^-1]`という省略記法を明確に訂正する。既存のsource・証明書・hashは変更しない。

**D_perpを(m_R0,w,positive odd)で完備化したw=0の環と読むなら、
その式をw=1のdehomogenizationに使うのは誤りである。**
例えばk[[w]]ではw-1が単位であり、w=1近傍の形式環はそこから得られない。

正確に使う環をD_algと書く。これはR0上の有限w多項式graded algebraであり、
CA側のdegree0 module C0、degree>=1 module w^d Hと、次数<=3の全actual関係から作る。
元B18原点のD_perpはD_algのw=0での完備化である。
一方、安定module M=H/Jinfについて成立する正確な式は

```
D_alg[w^(-1)] = M[w,w^(-1)].
```

この**未完備w多項式環を先に局所化**し、その後(m_R0,w-1)で完備化するとM[[w-1]]を得る。
nilmoduleと被約Nにも同じgrading操作を適用する。

D_algが単なる形式的な存在でないことも必要である。
実sourceのzero48単位Hensel系にはwが入らず、positive5線形系の係数と分母もzero sourceだけに
依存する。rank-one shearは正確に有限w多項式であり、全150出力はw次数<=3。
zero48の単位解はalgebraic power seriesなので、係数はpointed étale algebraic近傍で扱える。
同じunit分母がw=1でも使える。従ってこのsourceは実際のglobal normalization source chart上で
W1へ延びる。既受理の全非零Wにおけるproper normalization fiberの唯一の幾何点と、
同source chartの単位条件により、W1で別枝を落とすこともない。

旧W1の218式との結合は、この実際のw-polynomial sourceのW1 germに対して行う。
w=0の完備環そのものにw=1を代入しているのではない。
したがってtail bindingのliteral旧49/group17/cut20座標、全中心一致、単位minor、
nilcurveと同R0 mapの計算は影響を受けない。
W1のactual familyを介する上記の順序が、元の省略記法に代わる論理である。

依存するfinite-w sourceは
`AMBIENT_B18_ACTUAL_TOP_COEFFICIENT_PREFLIGHT_20261104_v63`、
`AMBIENT_B18_FULL_RELATIVE_INPUT_20261104_v63`。
W1 proper fiberと実交叉認識は既受理の
`AMBIENT_B18_ALL_NONZERO_W_FULL_COMPLETIONS_20261103_v63`と
`AMBIENT_W1_S61_INTERSECTION_RECOGNITION_20261103_v63`を使う。
新しいCAS、jet、点、環の同型の追加推測はない。
