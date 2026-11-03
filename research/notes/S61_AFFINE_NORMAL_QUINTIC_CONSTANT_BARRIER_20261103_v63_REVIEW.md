# Quinticのformal45定数cofactor障害：独立監査

2026-11-03。凍結 `S61_AFFINE_NORMAL_QUINTIC_CONSTANT_BARRIER` を**自由formal45係数環**の障害として承認する。
R'=K[n0,...,n41,S1,S2,H]、K=C(beta)とする。h次数0--4の係数を5成分vectorに取った全670式には、
parameterに依存しない非零K係数の行結合で一次以下となるものが存在しない。
非線形の実normal recipe S_l(n),H(n)へのbase change後に同じ負の結論が成り立つとは推論しない。

## 親quinticとの対応と別の670行構成

親packetの独立reviewをhash固定した。このreviewは元の全134式のcentering、実plane拘束、
全134のQQ代入・division恒等式、generic betaでの有限Macaulay cofactor recipeを検査済みである。
今回はその原式対応を引き継ぎ、全134 remainderを改めて全て読み直す。
元134式の一部だけを選ぶことや、線形化だけに置き換えることはしない。

作者の各h^j倍を再びlong divisionする方式に対し、独立実装ではmonic G=g(h)-Hのcompanion actionを使った。
vector=(v0,...,v4)へのh作用は

    w0=(H-g0)*v4,
    w_l=v_(l-1)-g_l*v4,  l=1,...,4

である。134個のvectorに4回作用して全670行を構成した。
全entryのformal parameter次数は2以下、定数項は零である。
全二次5145列と一次225列を別にし、QQの全分母が32003の単元であることを確認した。
全237811個の非零modular entryを元native入力と一つずつ照合した。

## 非零670次minorはrank下界の証明

作者の全670行・選択670列を、保存された(h次数,parameter monomial)のlabelまで一致させた。
独立C++ではFLINT rref/determinantを呼ばず、整数modular Gaussian eliminationとFermat逆元を使う。
全元係数から作った670×670行列式は正確に4745 mod32003である。
従って同じQQ minorは非零のp-unitであり、QQでの二次block rankは670。
上界670は行数なので等号である。modpのcolumn spanをQQへ持ち上げた論法ではない。

一般betaについては、親reviewで承認されたbeta有理plane frame・primitive h・有限cofactor逆行列を使う。
選択されたminorの係数はそのbeta-only開上の有理函数であり、凍結QQ中心で定義され、
非零4745というunit特殊化を持つ。従ってこの函数は恒等零でなく、一般betaでもrank670である。
未知qやnormal parameterの値に依存する局所化を追加していない。

## plane originに依存した偽の一般化を避ける

plane内でoriginをc=(c0,c1,c2)だけ動かしてh'=h-c0とすると、

    p_l'(h')=p_l(h'+c0)-c_l,  g'(h')=g(h'+c0),
    S_l'=S_l,  H'=H.

normal座標nは同じplaneを切るので変わらない。
元式のh^j倍から新しい(h-c0)^j倍への変換と、targetのpower basis変換は、K上のunipotent binomial行列である。
従ってformal parameter次数を保存し、二次係数rankを変えない。
独立コードは一般c0と一般g0,...,g4,Hを記号のまま使って、binomial逆積・det1と
新companionが C^(-1)(T-c0 I)C であることを全25entryで確認する。
generic originが既知五点の一つという根を選ぶ仮定は不要である。

## 結論の限定

二次blockの行核が零なので、parameter-constantの非零行結合は必ず二次部分を残す。
これが主張する純一次consequenceの零性である。hのcofactor次数0--4は全て含まれる。

S_l,Hは実代数ではnormal polynomial circuitsであり、自由parameterではない。
その非線形base changeはformal二次項のcancellationを起こし得るため、**actual42-normalだけの同じ障害は未証明**である。
parameter-dependent cofactors、より高いmodule closure、unit、一般五点尽くしも排除しない。
この監査では高次cofactor Macaulay、GB、探索の拡大を実行していない。

code: `code/s61_affine_normal_quintic_barrier_review_20261103_v63.{py,cpp}`。
receipt: `reports/S61_AFFINE_NORMAL_QUINTIC_CONSTANT_BARRIER_20261103_v63_REVIEW/`。
