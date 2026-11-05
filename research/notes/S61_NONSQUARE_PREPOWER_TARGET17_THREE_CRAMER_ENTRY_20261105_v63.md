# 三 Cramer 分子の worker custody

2026-11-05。SOURCE fcf73502 の同じ三分子だけを固定する。κ17独立peerとroot確認前に
polynomial workerを実行しない。det A3=-κ3、κ3逆=detA8/(cprod Saux T) の既存単元を
使うが、ここで逆元の係数を計算しない。

H5 certificate の全 JSON は読まず、必要な raw_compact14x5 と rho の subtreeだけを
選択的にdecodeする。六heavy入力のQQ定数scale、三anchor、三v17をliteralに固定する。
十二入力・九2×2cofactor・九term積・六accumulator加算からなる54 node DAGを使う。
cofactorのsignは行列内添字の(-1)^(i+j)だけで、Cramer列置換の再signは掛けない。

整数共通分母と暗黙格子添字を用い、各係数、sparse packed key、格子row offset、
denominator、poly headerとmetadataを全て計数する。上限73322、temporary434bit、
入力と新成果物合計5.82MB、60秒/2GBを保つ。全積は逆順再生、全加減算は全係数
backaddで検証し、三出力を一ファイルずつQQ既約化・再読込する。

得られるYに対してx前三=-Y/κ3、x4=0。普遍adj(A3)恒等式とκ17零から第四行も従う。
これは同じ列17のpair補正だけで、残る十角補正・最後のentry・rank・全FULLU57を
計算したとはしない。旧四＋十 sourceとtailは全て保持する。
