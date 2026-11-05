# lower四channelのdual行を直接構成する精密化gate

元の441bit priorでは、不要なfull B99逆行列と冗長なC*J検算が過大な正上界を作った。capは変えず、表現を変える。必要なのはs=h B99という四つの特定行だけである。normalized basis n_p=Ystar0 rep_pの先頭が1なので、free座標でs=0とし、p降順で `s_p=h rep_p−Σ_{j>p} n_pj s_j` と構成できる。これは新しいGaussianではなく、凍結されたunit基底の三角式である。

正のenvelopeはこの特定行の右辺、tail、back-add検算を全て含む。C*J=0は既受理F25 PiY0 J=Iから厳密に従うため、大きい別の検算をsource builderの要件にしない。F25の固定membershipは保持する。cap384を覆うかを先に確認し、覆わない場合は署名付き行やsource値を構成しない。

この精密化はsource係数の表現だけを変え、四つのactual lower formula、619のoriginal source、pivot carry、G3一次clearを変更しない。adjoint solve、X2値、新physical/gauge/jetは計算しない。
