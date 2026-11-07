# normalnil非可換affine型の全重複度上界：独立補足

2026-11-07、主義務G。別作者による全source幾何の確認。**PASS**。
`PLANE_RANK_ZERO_NILPOTENT_AFFINE_ALL_MULTIPLICITIES.md` の主張を受理する。
実Dが既約被約quintic、rank零planeのglobal-vector normalnil型でb非零なら、
**全重複度でambient projective像次元は22以下**である。

全m modular補題とその独立reviewはconormal可換性を位数二なしに与えるため、
全35/30/9源を使える。
先行独立review `PLANE_GLOBAL_NILPOTENT_NONZERO_B_SEPARATE_REVIEW.md` のJacobi、
normal shear/scalar、全centralizerと逆は位数二を使わず成立する。
したがってVはJordan3かJordan21、前者の全sourceは一つのR orbitでleading determinant非零、従ってm二である。

Jordan21では、全解は

`V=N, T=Dρ+xN+yB+vC, K=zw`,

`det(Ebase,T,N)=x1²*(v*x1−(1+3ρ)*x2)`。

従ってm≥3の全補集合は正確にρ=−1/3,v=0である。
このρでは2+3ρ=1なので、actual exp(sB)の非共鳴作用でyを消し、exp(qN)でxを消せる。
v=0は保存され、全逆は有限unipotentである。
よって全補集合はNDIAG(−1/3)という一つのwhole normal formへ入る。

元のrank一nil V orbit四、normalnil line一、Kのa,b二、残るT解x,y二を戻した
fixed-plane affine symbol baseの次元は九以下。
全source十一からρ固定とresonant row零の独立二条件を引いても同じ九である。
正規化点だけを数えたり、x,yをhigher fiberの自由度へ混同したりしていない。

先行の独立Fraction report
`CLASSIFICATION_PLANE_GLOBAL_NILPOTENT_NONZERO_B_REVIEW_20261107.json` は、このρ=−1/3も
当時のm二外として明記したまま、**全350×30/9**、nonzero maximal minors及びcomplete kernelsを既に検証している。
そのexact nullityは7,0である。
元producer reportの全350係数を新たに再実行する必要はなく、今回その同じ不変証明書と入力hashを確認した。

全higher解はPi1∈ker D1、任意Pi1上のPi2は同一D2のinhomogeneous fiberだから次元七以下。
全plane六、tensor射影化−一を戻せば

`image dimension≤9+7+6−1=21`。

m二の全branchは先行whole theoremで22以下なので、全mでも22以下となる。
旧certificateの `m2_excluded=true` をfalseへ書き換えず、新しい全m入口による適用だけを追加した。
normalnil b零、normalss高m、traceless symbol零、他co-Higgs型、全G/Eの尽くしは結論しない。
成分数を変えず、909原稿/PDF・global/Gitは変更していない。

全入力hashとpure-proof有限契約は `jobs/plane_nilpotent_affine_all_m_review.json`、
独立受理receiptは `reports/CLASSIFICATION_PLANE_NILPOTENT_AFFINE_ALL_M_REVIEW_20261107.json` に保存した。
新しいCAS、operator探索、jet又はGBは用いていない。
