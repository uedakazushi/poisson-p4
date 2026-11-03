# 全射影 CB の随伴 P1 模型と階数二境界：独立監査

2026年11月3日。frozen QQ producer の再実行と独立 sparse150 物理 tensor 検査を行った。
全22 raw/canonical 基底、全352 Lie 作用列、比率補正を含む304相対作用列、h-line を忘れる
二つの部分空間の追加 M01 作用、全22パラメータ raw Jacobi、rank1 minors と実 quintic を確認した。
三つの境界点の全像微分行列と非零 QQ minor27,28,24、全175×5定ベクトル因子行列も再構築した。

比率 µ=b/a を定数としたまま全 parabolic 作用がファイバーを保つと仮定していない。
実際の f/g 三係数 block は同じ binary 変換に対し、character
det(Dzw)/(cx cy) だけ異なる。独立 Lie 作用で差χI3と χµ∂µ の補正を確認した。
連結 parabolic 上ではこの有限 character が Lie 差を積分し、非自明な随伴 P1 束が定まる。
その上の tautological 条件 (f,g)∈(a,b)⊗Sym²W と H16 が rank19 の閉部分束を与える。
従って滑らかな projective source の次元は旗9+比率1+projective fiber18=28。
全 raw0i=0 から Euler-vertical kernel は零であり、projective map は basepoint を持たない。
a,b≠0では character が全 Gm を動かすので比率を1へ移せ、既存 CB の dense 族と一致する。
properな既約閉像は全被約 CB となる。

raw Jacobi の唯一の成分は g f_x−f g_x で、その三係数は2×3 rank1 idealの三 minorそのもの。
従って rank1 cone の [a:b],F による解消は全 Poisson source を与える。
実 quintic a b h F² から rank≤2 逆像は a=0,b=0,F=0。
F=0像は比率を忘れて両 endpoint に入るので、境界は二つの proper既約像で尽くされる。

a=0 と H-only 空間は旗の H-line を忘れる17次元連結 parabolic に保存される。
元の16作用と追加 M01 に対する不変性から、底は dim8の旗(L,W)となり、projective像上界は26と23。
b=0は元の旗dim9を使って上界27。
独立に再構築した affine 微分minor27,24,28がそれぞれの下界と一致する。
b=0全像には定ベクトル因子がある一方、a=0保存点の intrinsic因子行列はrank5である。
したがって二つの像は相異なり、被約 boundary の全既約成分は正確に26次元と27次元の二つである。

これは generic Poisson rank≤2 の境界を含む全射影 CB の支持と次元の結果である。
全正規化、全 R1、CAとの交叉、各境界点の旗ファイバーを結論しない。
code/ambient_CB_global_flag_review_20261103_v63.py と同名 REVIEW report が独立証拠を保存する。
