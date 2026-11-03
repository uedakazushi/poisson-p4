# 受理済み実交叉定理からの帰結：埋込み associated prime の存在

2026年11月3日。対象は S61 内の D=CA∩CB であり、個々の CA,CB,N は被約とする。
既に独立監査した次の三事実だけを用いる。

1. Dred=N は既約である。
2. rank Kβ=3 の非空開部分では D は滑らかな Z3 上の affine9 bundle であり、被約である。
3. 指定 a=1 Q および共通10次元非正規層の全点で D は非被約である。

従って D は generically reduced だが、その各非被約点 x における Noetherian 局所環 A=O_D,x は
少なくとも一つの埋込み associated prime を持つ。

実際、P=sqrt(0) は A の唯一の minimal prime である。仮に埋込み associated prime がなければ
Ass(A)={P} となる。P 外の各元は非零因子なので A→A_P は単射である。
一方 generic reducedness により A_P は体であり、A の全 nilpotent はその中で零となる。
単射性から A は被約となり、第3項に反する。

これは単なる正規化 tensor product の embedded structure についてではなく、実際の target
scheme intersection についての結論である。ただし埋込み prime の具体式、個数、次元、
nilpotence order、全 primary decomposition は決めていない。
元の射影空間についての帰結も指定点を含む generic Poisson rank4 開部分内に限る。

入力は AMBIENT_FULL_FIRSTJET_FIBERS_20261103_v63_REVIEW.md と
AMBIENT_INTERSECTION_20261103_v63_REVIEW.md、および同名の独立検査 report である。
追加計算を要しない Noetherian 局所環の帰結として記録する。
