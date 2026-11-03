# Controlled Schur unitの代わりの特定cubic-power membership

2026-11-03。一つのgeneric aとsharp532 pivot delta0をdeltaと略す。670本のSchur scalar generatorsのidealをJ_gとする。未知残余が空であることは、単なる不特定のlocalization powerではなく、**特定のmembership delta³∈J_g**と同値である。membership自体はまだ証明していない。

## 1. universal determinant identity

Bを任意の5×4 pivot matrix、Xを任意の5列とし、delta=det B_bottom4、g_j=det[B,X_j]とする。top rowを

    delta*top−B_top*adj(B_bottom4)*bottom

へ変えるとXのdeterminantはdelta倍となる。変えたtopのj成分はg_jである。Laplace expansionから

    delta*det X=sum_j(-1)^j*g_j*det(X_bottom4_except_j)

を得る。producerはBの20変数、Xの25変数を全て独立なinteger variablesとして、両辺の全2880非零係数を正確に検査する。

XにPhiの任意の五列を入れれば、実idealとして

    delta*I5(Phi) subset J_g

となる。pivot columnsが元Phiの定数結合でも成り立つ。未知supportについての近似ではない。

## 2. 必要性

Bが空ならactual Aはknown E5だけで、formal R'上のmoduleは(R'/m)^5である。従ってFitt0(A)=I5(Phi)=m^5となる。全Phi entriesはmに入るのでdelta∈m^4、delta²∈m^8⊂m^5である。第1節を使うと

    delta³ in J_g

が従う。

## 3. 十分性

J_g⊂I5(Phi)⊂Ann(A)である。delta³∈J_gならdelta³はB(a)をannihilateする。しかしsharp moment theoremでdeltaはB(a)の非零因子であり、B→B(a)_deltaはfaithfully flatである。従ってB=0となる。

以上からexactに

    B=0 iff delta³∈J_g

である。known-five factorはdeltaで自動的に消え、unknown associated/nilpotent partsはこのcontrolled localizationで失われない。

## 4. 範囲

deltaのformal parameter次数≤8なのでtarget delta³の次数は≤24である。ただしinhomogeneous generatorsを使うmembershipのcofactor次数には、ここから上界を与えていない。degree≤14のcofactorで足りる等とは推論しない。新しいMacaulay/GB探索、巨大なactual deltaの展開、actual membershipの確認は行っていない。

これはgeneric five-point exhaustionの明確な有限algebraic targetを与える補題であり、尽くしそのものの証明ではない。`code/s61_affine_normal_schur_cubic_power_20261103_v63.py` と `reports/S61_AFFINE_NORMAL_SCHUR_CUBIC_POWER_20261103_v63/` に、author/independent source hashesとuniversal identityを保存する。
