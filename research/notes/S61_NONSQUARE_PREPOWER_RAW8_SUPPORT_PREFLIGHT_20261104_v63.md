# 同じ compact7・target19 を raw8 minor に戻す source 前見積

2026-11-04。scalar値、minor、rank、点、primeは新たに計算しない。
候補は以前のborder6で固定した compact7/元target19だけである。

元prePOWER15×22のrow順 `[0,1,3,6,2,4,5,8]`、元32column順
`[22,11,13,25,15,19,27,31]` の8×8をA8とする。先頭四行はpureで、後の四行は
compact `[1,3,4,7]` のpair-only sourcesである。全sourceとq-DAGの実rootを保持する。

二つのpair-kernel行の2×2cofactorについてJacobi complementary-minor identityを
用いると、先行cleared border6 determinantは

    detborder6 = rhoLCM4*rhoLCM7*delta2 /
                (rho1*rho3*rho4*rho7*gamma) * det(A8)

となる。両側を既知FULLU unitsだけで分母消去すれば多項式恒等式であり、最初に
D(delta2)で導いてintegral domainへ戻せる。従ってこの恒等式自体はpを反転しない。
det(A8)を直接扱うことで、delta2での未知除法や7022項Sauxの掛算を必要としない。

さらに既存B3 rows current15 `[7,13,10]`、heavy targets `[28,30,29]` を加えると、
先の八行のこの三列はliteral zeroである。11×11はblock lower triangularになり、
determinantは `det(A8)*Delta407*(known row rho units)`。もしこの固定det(A8)が
非零と後で証明されればD(det(A8)*Delta407)上でsame-cokernel4×11へ進める。
残るcurrent15 `[9,11,12,14]`、pure targets `[1,2,3,4,6,7,8,9,14,17,18]` で、
全8Q9tailを保持する。ここではdet(A8)の値も非零性も主張しない。

計算recipeはtop4行を選ぶpure6列の15通りLaplace展開である。各項のbottom4 determinantは
共通のheavy2 minor六個とpure2 minorから作る。top4 minorは一つずつ計算・消費し、
全15個を同時保持しない。元P4の一つだけは受理済み818項D4をreuseする。各rowの
全入力monomialのa指数下界を取り、そのknown a冪だけを予め除く。Cやpの未確認因子を
仮定して除かない。全atomのa/QQ denominator形状を検査し、B poleが現れたら停止する。

totaldegree interval＋exponent boxの整数格子点を全nodeで数え、l1=0以外のcancellationを
仮定しない。last-use解放、全literal atoms、output-sized scratch、15000reserveを含める。
QQ係数heightと分母、serialization上界も保存する。guardを通った場合でも、これは
source前見積であって新しいscalar計算の実行許可ではない。
