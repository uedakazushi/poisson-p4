# T37 low-degree leaf theorem：simple component の正確な specialization scope

2026-11-08、root手読解に基づく明示的射程限定。元 note `T37_WHOLE_CLOSED_LOW_DEGREE_LEAF_SURFACE_THEOREM.md` SHA256 `5a95d021d4a0a1ce998c61192fba47ad65605273b38cf99fe8d713e5ac6a40a7` と元12-input manifest の bytes は変更しない。

## 採用する必要条件

rank四 tensor π0 が whole T37 に属し、actual quintic D0 の **multiplicity一の prime component** B 上で、他の全componentと B の singular locus を避けた point p0 を取る。さらに π0(p0) の rank が二とする。このとき p0 は degree≤5 の integral Poisson invariant surface に属する。

元 note 第3節は D0 自身が p0 で smooth という仮定を明記している。第4節の「各 component の一般 smooth rank二点」は、上記の simple component 上の D0 自身が smooth な点に限定して用いる。multiplicity二以上の component の reduced support だけが smooth な点へ、本証明をそのまま適用する主張はしない。

## generic point の選択

curve-selection の DVR family で D_R は p0 近傍において smooth of relative dimension三。完備化した局所環は R[[u1,u2,u3]]。generic の pencil base／bad level／singular locus は proper closed subset であり、その非零局所方程式の pullback はこの formal ring で非零である。u_i=uniformizer*α_i とし、有限個の最低重み係数 polynomial を消さない一般 α_i∈C を選べば、p0 に特殊化しつつ generic bad locus を避ける section p_R を得る。必要な有限 field extension は元証明どおり許す。

これで generic leaf surface Z_K の Hilbert properness、flat ideal の Hamiltonian stability、char零 radical/minimalprime stability と rank二による低次元 component 排除の全議論を、元 note の第3節の仮定内だけで使える。

## 現在の target への適用

Ga9 の actual D0=−4hF4 は被約、h≠0 上で F4 は w 一次の smooth hypersurface。第2節の全 leaf 計算は h=1 にあり、ここで rank二。従って `GA9_LINEAR_W_WHOLE_T37_LEAF_DEGREE_EXCLUSION.md` の degree九対五の矛盾は限定後も全く変わらない。

同様の simple cubic target へは、それぞれの単純 prime component の ambient-smooth rank二 open 上でのみ使う。nonreduced component 上の一般点被覆、normalization の flatness 又は新しい numerical computation は追加しない。
