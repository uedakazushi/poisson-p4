# 平行 Lagrangian jet の普遍特異点と六成分排除の独立監査

2026-11-03。`Q_LAGRANGIAN_JET_SINGULARITY_20261103_v63` の結論を承認する。
指定した六つの全 actual 成分の一次 jet discriminant は恒等的に零となり、
正 weight6 軸の全非零点はそれらの外にある。MPS 成分一般や未知成分への拡張はしない。

まず仮定された幾何形を持つ線形 jet の全12係数から出発した。
constant J との compatibility の全四式は rank2 で、fiber-fiber の二係数を一意に消去して
全10係数の compatible 空間を得る。primitive perspective 補正後の24係数の像は rank8 であり、
表示された八変数族の係数空間と、結合行列の rank まで含めて厳密に一致する。
したがって特殊な ansatz だけの性質を全 jet に広げたものではない。

八変数を残したまま primitive beta、標準 Pluecker 行列、canonical covector、5×5 skew 行列を再構成した。
全五つの符号付き Pfaffian が (Pl*covector,Pf(Pl)) に一致し、
点(0,0,0,0,1)で全式が零、全 Jacobian は一列だけが非零で、その最終成分は定数−1となる。
ambient Pluecker quadric 自体はこの点で滑らかだが、その上で曲線の残る全微分が零であり、
曲線の接空間次元は少なくとも3である。従って degree5 Pfaffian genus-one model は滑らかではない。
既存の Fisher–Sadek discriminant criterion により Delta=0 が恒等的に従う。
特殊化でモデルの次元が上がる場合も、Delta≠0なら滑らかな genus-one curve という定理に反するので同じ結論となる。

rank4 の任意の値 [[0,A],[-A^t,B]] も全五変数のまま正規化した。
表示された block lower-triangular frame の M*C*M^t は標準 J、detM は1/detAである。
vertical plane を保ち base 座標を変えないため、仮定の形は値の正規化で失われない。
平行移動もこの形を保つ。primitive 補正は base の線形式だけを使い、vertical なベクトルを加えるので、
普遍八変数 jet の適用に未検査の rank4 点は残らない。

六つの凍結全 symbolic source について、実際の projective affine bivector
Pij−xi*P0j+xj*P0i を別に作り、P_xz=0 と全八つの fiber 偏微分の零性を全係数で確認した。
特に P0i を持つ C+ source の単なる4×4部分行列を使用していない。
既存の全 source の PGL 稠密性と、primitive first-jet curve の frame・scale に依存しない特異性を使えば、
Delta の零性は各成分の全 S61 閉包へ延長する。単一標本の消滅から普遍恒等式を推測していない。

正軸 eta1 の beta は実 S61 再構成と affine 微分から再計算した。
既存の厳密 Fisher–Sadek 演算で c4=0、c6=−949978046398464、
Delta=−522255954073519803473068032 と確認し、Delta=(c4³−c6²)/1728も照合した。
非零性は六成分の普遍零性と矛盾する。既承認の有限 diagonal 作用により全非零 eta は一軌道なので、
この排除は軸の全非零点に及ぶ。

fresh producer と全証明書も一致した。独立コードは
`code/q_lagrangian_jet_singularity_review_20261103_v63.py`、receipt は
`reports/Q_LAGRANGIAN_JET_SINGULARITY_20261103_v63_REVIEW.json`。
