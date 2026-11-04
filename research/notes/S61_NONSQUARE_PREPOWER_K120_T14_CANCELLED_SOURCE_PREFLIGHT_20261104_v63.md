# gを多重線形性で先に除くK120/target14のsource費用

2026-11-04。新しい係数計算を行わない。選択七pure列
[22,11,13,25,15,19,14]のraw pair比例は、κ14だけでなく三raw交叉の零で既に確認済みである。

A=current11/c10、B=(current12−current10)/c11、C=current7/c6 とし、
K115=gA−hC、K120=gB−h120Cを使う。最初五行P4,G20の下へこの二行を追加したD7は

    D7 = g² det(P4,G20,A,B)
         −g h120 det(P4,G20,A,C)−g h det(P4,G20,C,B)

となる。Cが二度現れる項は交代性により恒等零である。
従ってgでの数値的除法を行わず、三項を保持した分子N7=D7/gの有限polynomial recipeを得る。
残余G120/target14はc11 N7/(a² N6)。分母のaとN6は現在のchartで既に単元であり、
g逆・h22逆・p逆・R434逆を追加しない。tail2は以前のpivotと全逆変換に残る。

費用は上五行の五minorと下二行のweighted二minorに沿う二十一項で見積もる。
後者はg(AjBk−AkBj)−h120(AjCk−AkCj)−h(CjBk−CkBj)であり、三項を削除しない。
各項のQQ分母と整数l1上界、degree/box格子、sourceの最後の使用、scratch、
元atomと最終整数格子の同時liveを含める。最終二十一積だけを一つの整数格子へ融合する
recipeについて、全scaleの整数性を検査する。実値・非零・rankは未計算である。

上限を超える場合は、その表現の失敗を保存し、係数計算やcap増加を行わない。
raw十minorの以前の上限超過記録も変更しない。全補集合と残る十targetを保持する。
