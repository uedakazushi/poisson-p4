# NEW firstjetのgeneric Sp4 stabilizer：独立centralizer監査

2026-11-04。対象は frozen `Q_NEW_FAMILY_SP4_CENTRALIZER_FREE_20261104_v63`。
使用するのは同じ旧OCTIC_BASE_DOMINANCE index0のQQ βとprime32003だけ。
元physical quadratic covariantsの全係数・同時共役性と実source点への束縛は、
既受理のSpin5/old-pointレビューをhashで固定する。新しい点やprime、rank探索はない。

C,Dの計4352 quadratic係数を全て読み、その分母と旧βの分母がp-unitである
ことを確認する。全32値をQQで評価し、係数とβを先に有限体へ写した全32値と
独立に比較する。保存されたC,DのQQ値・modp値・旧Krylov値とも一致する。

32×15のcommutator mapは作者の行列積routineを使わず、matrix unitの式
`[E_pq,C]_ij=δ_ip C_qj−C_ip δ_qj` とdiagonal差の式から別に組み立てる。
全480成分を保存値と比較する。使う15行は既存certificateで指定されたものだけ。
その15×15 residue matrixのdeterminantを、作者のmodular eliminationではなく
整数lift上のfraction-free Bareissで計算し、mod32003で **17760≠0** を確認する。
別のminorは選ばない。

各entryはβのp-integral homogeneous quadratic polynomialである。従ってこの
指定minorはdegree30の実QQ polynomialで、その **同じ実QQ点**での値がp-unit。
ゆえにQQ上でも非零、NEW像B上で恒等零ではない。modp rankから未制御の
char0 fiberを推定する議論ではない。非空なrank15 openが得られる。

このopenの任意の幾何点で、C,Dと可換な4×4行列Tについて
`T−tr(T)I/4`はsl4のkernelとなり零。よってLie代数だけでなく **全matrix
commutant** がscalarである。Sp4のβ-stabilizerは共変性によりこのcommutantへ
入り、symplectic scalarは±Iに限られる。元W16は三つのcovariant indexを持ち、
中心−Iはβを−βへ写す。rank15ならβ≠0なので−Iも除かれ、stabilizerはidentity。
線形kernelの議論は代数閉包への拡大後も有効であり、char0のstabilizer scheme
についても自明となる。有限・非連結stabilizerをLie計算だけで落としてはいない。

選んだ一つのminor chartをSp4安定とは仮定しない。全rank15条件は同時共役に
不変で、その非空openをNEWの既約性と合わせてgeneric freenessを結論する。
これは **firstjet B_NEWのgeneric作用だけ**。全点のfreeness、25次元Poisson
成分、target field index、orbit section、B_NEW自体の有理性を本packetだけで
結論しない。元wrapperの算術前repair履歴も保存されたまま照合する。
