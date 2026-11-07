# nil31短鎖 r=2/3：実三超平面・rank四quadricの二flag診断

2026年11月7日。E_multiple_hyperplane。全源の帰属は未解決。
これは既知MPS閉和への二つの特定の実flag族を調べた負の診断であり、
任意のMPS変形や成分帰属を否定する結果ではない。

既済の全六変数源は `MODULAR_NIL31_SHORT_SUPPORT.md` §7 により
二変数 t,d の模型へ正確に移り、td≠0は一つの22次元PGL軌道となる。
従って既知閉和の中に t=d=1 への実曲線を一本与えれば全源を被覆できる。
この最弱目標を保持した。

入力親族は `THREE_HYPERPLANES_QUADRIC_INDEPENDENT.md` §2 の
独立三超平面 hs,X,Y と Q=XZ+YW の全源である。

    D1=X∂X+W∂W, D2=Y∂Y+Z∂Z, N=Y∂Z−X∂W,
    P=a D1∧D2+b Q∂Z∧∂W+U N∧D2+V D1∧N.

a,b,U,V は本当に独立な四係数である。a≠0なら元のp=U/a,q=V/aに戻り、
N∧N=0なので係数の積を独立変数として誤って自由化したものではない。
任意四係数のJacobiは多項式恒等式であり、既知MPS閉和への帰属も稠密ab≠0から閉性で延びる。

第一flagは

    (hs,X,Y,Z,W)=(h+εy,h,z,y+εx,−w).

行列式は非零定数倍のε²、逆行列のpoleは二、bivectorのpoleは四である。
四実係数を各ε次数四以下の多項式として全20係数を保ち、
Euler鉛直を除いた全150係数で負冪の消失と端点を比較した。
系は63×20、rank10、拡大rank11。
保存した二行functionalが行列を消して右辺2/15を与える。
実quadricのvertexは e_x−εe_y で、target cofactorの長さ三鎖の次段を動かさない。

第二flagはこの幾何を変えて

    (hs,X,Y,Z,W)=(h+εy,h,z+βεy,y+εx,−w)

とした。vertexは e_x−εe_y+βε²e_z となる。
βを標本値で掃引せず、QQ(β)上の全69×20系を厳密消去した。
一般rank10、拡大rank11で、非零残式の分子のgcdは1。
選択minorの例外はβ=0だけで、この点も別の有理functionalにより不整合となった。
従ってこの一つのβ旗族にも、指定した四実係数の正則多項式族による端点は存在しない。

exact producer/report/jobはそれぞれ `classification_nil31_short_flag_confluence` と
`classification_nil31_short_chain_flag` のstem。実行0.318秒と0.373秒。
全行列、右辺、非零minor、例外処理、実行sourceと入力hashを公開reportに保持した。
gate PASSは契約の整合確認であり数学の独立レビューではない。

両方とも同じ次数を増やして再試行しない。この二flag経路を停止し、
既知25次元成分の真正graph認識又は別の全源構造へ切り替える。
任意Laurent親係数、異なるflag、任意MPS族の不存在は主張しない。

## 別構造による二つの有限確認

同じ稠密軌道のt=d=1を、実symplectic点 `(h,x,y,z,w)=(1,0,1,1,0)` で
Darboux正規化し、忠実S61の全150係数再構成・全134元方程式を確認した。
T39とT40の普遍被約graph全45式は、双方とも成立しなかった。
`CLASSIFICATION_NIL31_SHORT_ELLIPTIC_GRAPH_20261107.json` に実frame・全61座標・
全90残式を保存した（3.976秒）。別の評価点へ取り替えて再試行しない。
これはこの既知二成分の除外であり、MPS40全体の非所属ではない。

またquadric Q=hy−zw自体を平滑化する、独立した構造も一度試した。

    Qε=hy−zw+εx²,
    Xε=(−2εx, y−15h, z+30εx, 0, h).

厳密に Xε(Qε)=0、`Xε⁵=60ε Xε³` である。
ε≠0で `S=Xε³/(60ε), N=Xε−S` は可換Jordan部分、SN=0,N³=0であり、
rank二半単純場とnil3からなる既済smooth-quadric mixed源へ属する。

この源の移送をdivisorの形だけで判定しなかった。
Q=EᵀBE、U,VをSの零・非零固有部分の射影、EU=UE、EV=VE、
qU=Q(EU)、qV=Q(EV)、KN=NB⁻¹、KS=SB⁻¹ とすると、全一period源は

    −5qU KN/2−2 EU∧NE+EV∧NE/2−EU∧SE/3
      −5qU KS/6−qV KS/2 + m SE∧NE.

この公式は公開mixed完全源のm7=5/6と全150係数で一致する。
今回の平滑化ではparticular項のε付値−1、period tensorの付値0なので、
正則端点に影響する自由係数は `m=c−1/ε+c0` の二つで尽くされる。
全55係数式はrank2、拡大rank3で不整合。
一行functionalが右辺−5を与え、periodで変更できないh²係数が原因である。
`CLASSIFICATION_NIL31_SHORT_SMOOTH_QUADRIC_BRIDGE_20261107.json` は
全公式・Jordan場・行列・右辺・実行bytesを保存する（0.590秒）。

この指定されたquadric/cofactor平滑化だけの負の結論であり、
任意Slodowy方向や一般平滑化を排除しない。
未完了Eはこのまま保持し、追加のcofactor/profile探索へ進まず、
root割当のrank零plane全支持の独立レビューへ切り替える。
