# κ14とraw比例を区別するsource gate

2026-11-04。current15のpair四行[2,4,5,8]を使う。新しい列14の係数はまだ計算しない。
既存κ四係数はexact公開値と独立peerへ固定し、元sourceの列14だけを必要部分へ制限する。

A=[Hbar27,Hbar31,M22]のsigned三cofactorがκなので、κ14=κ M14=0は
M14がAの三列の張る空間へ入ることを示す。M14とM22のraw比例とは異なる。
旧D(Saux T)上ではκ3が既に単元だから、alpha,beta,gammaを最初三行のadjugateで
有限に定義できる。κ14=0なら全十五source行に対し

    column14' = column14−alpha Hbar27−beta Hbar31

という可逆列操作を行えばpair四vectorはgamma M22となる。その場合でも、K115/K120の
heavy成分による補正は保持し、未補正の裸七borderへ置き換えてはいけない。
κ14非零は追加pure方向を示すが、既にextra2をpivotにしているため新しいrankではない。

raw比例を直接証明する十分条件は、既知非零anchor M[2,22]を用いた三式

    M[2,22] M[i,14]−M[i,22] M[2,14]=0, i=4,5,8

である。これらが全て零ならfraction fieldでraw比例し、全交叉minor恒等式は整域へ降りる。
一つでも非零ならraw比例は偽だが、κ14の零性とは別である。

この論理差を明示し、κ一式・raw三式・両者の三つの候補についてsource入力、
QQ高さ、最後の使用、scratch、元atomのreserveと保存byte上界を別々に記録する。
ここではどのscalar計算も実行しない。全tail2効果と全七remaining tail、
既存K120 raw十minorの上限超過記録を変更せず保持する。
