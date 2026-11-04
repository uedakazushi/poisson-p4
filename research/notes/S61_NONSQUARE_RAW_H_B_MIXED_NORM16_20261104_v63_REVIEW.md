# mixed m1 norm16の独立レビュー

作者のscalar Bareissコード・builderはimportせず、保存sourceの全3 Mi=Gmi、
選択m1のQQ normalization、Hのeモニック16・全16 source除算を独自整数演算で検査する。
同じ256entryをQQ定数で整数化した行列について、別のFLINT polynomial-matrix
library determinantを一回だけ呼ぶ。solve・他norm検索・点・素数・GBは使わない。

native determinant全係数とdegree98/a-order35を作者出力と比較し、QQ content・
a35除去・g13の両quotient・巨大整数Bezoutを全係数で検査する。g13は必要a条件で、
非定数であることは実際のbad存在を意味しない。

同じh1基底なので、R=QQ[a,a⁻¹]上T=R[e]/(h1,k1,m1)は横結合[Mk1|Mm1]のcokernel。
Fitt0(T)は両determinantを含み、既知単元を外したgcd g13も含む。PID上のSmith積を
sとすれば(g13)⊂(s)、従ってsはg13を割る。Tは旧有限QQ代数の商で、QQ-length≤13。
実際のdモニックH4を加えるcoverはlength≤52。

このreviewのactual Bbadへの適用はD(G)側。G=0のQ/Lを他のallpoint unit証明で
除去しない限り保持する。全域emptyやfull cubic/quartic結論は追加しない。
