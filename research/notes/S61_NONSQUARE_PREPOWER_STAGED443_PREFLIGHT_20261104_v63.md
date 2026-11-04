# 同じ11pivotを4＋4＋3へ分けるsource adapterとG115一行の費用

2026-11-04。新scalar parameter係数、minor、rankは計算しない。
source順はpair-only四行 `[2,4,5,8]`、pure四行 `[0,1,3,6]`、B三行 `[7,13,10]`。
target順は `[27,31,22,2]`、`[11,13,25,15]`、`[28,30,29]`。列単位はκ11と同じ。

A4=[Hbar0,Hbar1,M22,BM2]のsigned cofactor定義からdetA4=−T。
Adj4(A4)の最後の行は−κであり、その四entryは既に計算済みである。
残る12entryを展開せず、inverse=−Adj4(A4)/Tをlazy finite operatorとして保存する。

元8minorからpair-firstへ動かすrow permutationの交換数16、column permutationの交換数16は
ともに偶数。従ってdetA8=detA4*detQ4であり、κ11で証明済みκ3D8=cprod*Saux*Tから

    detQ4=−cprod*Saux/κ3,
    Q4^-1=−κ3*Adj4(Q4)/(cprod*Saux).

を得る。κ3はD(Saux*T)で既に単元と証明済みなので追加chartではない。

さらに、既知六pure列rank一により、u=P4[:,22]、A=P4[:,11,13,25,15] とすれば
Q4=A−u*v^t/κ3。v_jは、anchorをjへ替えた同じsigned三cofactorのκ3成分である。
tildeQ=κ3*A−u*v^t の全3minorには普遍的κ3²因子があり、そのquotientは

    ±(κ3*det(A_minor)−v_minor^t*Adj(A_minor)*u_minor)

というpolynomial recipeである。従ってQ4 inverseは

    −AdjReduced/(cprod*Saux)

となり、κ3を反転するnodeも新しいexact polynomial divisionも不要になる。
このrank-one cofactor identityは別の形式整数多項式演算で全16entryを確認する。
二つの4inverseと既存B3inverseを合成し、元15×22とのsame-cokernel4×11と両側逆を保存する。

## 一つのG115 source行

残る四B-coupled行のうちcurrent15=11、compact14=10、actual135=115を選ぶ。
既存B3 Cramer座標は `[lambda415,0,0]` で、対応pivot sourceはcompact6/current7/G85。
両rhoはc_i*L、L=aFC²なので、B三列を殺すprimitive sourceは

    K=Delta407*G115/c10−lambda415*G85/c6.

この二sourceだけを使い、三つの大きなadjugate座標を作らない。target19は既知六pure列の
一つなので、pair eliminationによるそのpure quotientはP5 pure blockを使う商と同じ。
従って実remainingrowのtarget19成分はc10/(Delta407*Saux)倍の

    det6(P4,Bhat4,K ; columns=[22,11,13,25,15,19])

としてsource-boundに指定でき、Tや新tail2の巨大積を使わない。値と非零性は未計算。

さらにsource-onlyの小さな別表示を保存する。g=N6,1、h=N10,1 は既存26/30項polynomial。
既存Cramer identityとdomain cancellationからg*N10,j−h*N6,j=0なので、
Ksmall=g*G115/c10−h*G85/c6 もB三列を殺し、Delta407*Ksmall=g*Kである。
g≠0の追加actual coefficient chartでは、同じtarget19成分をc10/(g*Saux)倍の
det6(P4,Bhat4,Ksmall)で指定できる。g=0を捨てず、三middle係数chartがD(Delta407)を
覆うことも保持する。一般field unitと全FULLU unitを混同しない。

二つのdet6は同じ一source行・一targetのrepresentation比較であり、支持・degree・source
profileだけを読む。whole-worker height/lifetimeが別に通るまで実計算しない。
全八tailのうち2は単位pivotで消費済み、残る七tailは4×11とfull inverseへ保持する。
実full57 source liftは別の凍結operatorに結び、projected sufficient criterionのscopeを守る。
