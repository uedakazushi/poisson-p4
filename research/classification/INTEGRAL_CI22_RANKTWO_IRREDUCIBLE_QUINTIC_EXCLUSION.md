# 全integral CI22・rank二導手曲面のactual既約quintic枝は空

2026-11-08、PRIMARY E。既受理全sourceと新しい有限因子証明を結ぶ純粋なaggregate。
新計算、source probe、GB、minor又はjetは無い。

**定理。** πをP4上の一般rank四Poisson構造、Dをその実被約既約quintic退化因子とする。
非正規曲面支持のintegral成分 S⊂SingDがCI(2,2)で、π|Sの一般rankが二となる場合は存在しない。
Sのnormality、零divisor Rの被約性、ordinary node、clean、vertical無し、normalization-flatを仮定しない。

## 全sourceの被覆

Sは非退化非planeなので、secant argumentによりDのS沿いgeneric multiplicityは二。
rank二のtangent blockが可逆でmixed blockがI_Sに入ることと、actualPf∈I_S²から
intrinsic conormal bracketは可換になる。
従って既存のfull150 tangency／intrinsic conormal全源が適用出来る。
Euler24は正確に商し、全free-coordinate inverseが保存されている。

Sがnormalなら `NORMAL_CI22_RANK_TWO_SEPARATE_REVIEW.md` の全空性定理を使う。
nonnormalなら、任意の複素symmetric pencilの全congruence分類とsize五・integralityの全補集合から、
`NONNORMAL_CI22_RANK_TWO_SEPARATE_REVIEW.md` §2の次の五actual surface orbitだけとなる。
generic pencilの分類を全pencilに外挿していない。

| 全surface model | 全Euler商source | 全actualPoisson支持での結論 |
|---|---:|---|
| regular22_1 | 10 | 二wholemapと二つの全五源zero補集合。F=0又はs²x4 q2² |
| nil32 | 12 | 全九inversechart。六mapはF=0、二mapはH+2Q、残A開はsmoothQ*C3 |
| singularM1 | 10 | 元350Jacobi ideal内でF=x0*(aq1²+bq1q2+cq2²) |
| cone_node | 14 | 全actualD係数が元Jacobi行のQQ線形結合、F=0 |
| cone_cusp | 16 | 全二五源とboth-zero補集合。decomposable又はprojective有理Casimir、F=0 |

表のsourceは必要supersourceを含む。残るJacobiを無断で零又は自由と仮定せず、
各actualPoisson点について以下の全源証明を用いる。

## 五modelの正確な結論

regular22_1は `CI22_REGULAR22_RANKTWO_WHOLE_FACTOR.md`。
旧独立inversecoverのt≠0全二parameter、t=0/p1≠0全三parameter、残二全五sourceを保持する。
rank四点のactualDはH+2rank3Qなので、指定する実既約被約Dでは空。

nil32は `CI22_NIL32_RANKTWO_WHOLE_COVERAGE.md`。
ここはq1=2x0x2+x1²+2x3x4、q2=2x1x2+x4²の**rank二full12 source**である。
第一mapの唯一Kを、既に単位であるA=p10だけで解き、残八mapと合わせた。
A開ではQは常にsmooth、C3の非零性も保たれ、actualDは全点でproper Q*C3。
generic C3の既約性と既約Gm×A3 domainの閉性により、この全mapは既知MPS閉和にも入る。
残二rank四mapはH+2Q、他はrank二以下。
rank零six-sourceの座標又は旧部分源を全十二源の代わりに使わない。

singularM1は `CI22_SINGULAR_M1_RANKTWO_WHOLE_FACTOR.md`。
全126actualPfaffian係数の元350へのQQ liftから、source係数を一つも反転せず超平面x0因子を得る。
全support場合分けを解いたこと、distinct-Q境界の全点membershipまで得たことは結論しない。
actual既約Dの空性にはその全H因子だけで十分。

cone_nodeは先行別作者review §7の全33original-row functionals。
cone_cuspは `CI22_CONE_CUSP_RANKTWO_WHOLE_EXCLUSION.md`。
全二openに元350と全126を戻し、both-zero全sourceも保持する。
いずれも全rank四点は存在しない。

全場合で、F=0ならrank四に反し、Fが非零のproper積なら実既約被約quinticに反する。
これで全normal／nonnormal integralCI22を尽くし、定理が従う。
pointwiseなfactor／空性の証明であり、Gの像次元≤24だけをEへ読み替えた結論ではない。

## 分類への接続と境界

全PGL normal formを戻してもrank、actual因子型、Poisson条件は保たれる。
零parameter、全inversechart、extra conormal section及び全Rを含むので、
未知の特殊β上の点をgenericなsource値の仮定で除いていない。
以上により、この全CI22 rank二conductor枝をEから除ける。

rank零conductorは別に既受理
`E_DEGREE_FOUR_RANKZERO_IRREDUCIBLE_QUINTIC_ALL_EMPTY` を用いる。
両者を合わせると、restriction rank零／二のいずれにもintegralCI22 conductorは残らない。
rank零nil32のsix-sourceと本稿のrank二twelve-sourceは異なる入力である。

projected VeroneseはCI22でなく、cubic scroll、quadric cone、他のfactor型、全E又は総成分数の
完了は本定理から主張しない。新componentの増分も零。
本aggregateの全依存SHAは同名manifestへ固定し、既存source・証明書・主稿・global台帳・Gitは変更しない。
