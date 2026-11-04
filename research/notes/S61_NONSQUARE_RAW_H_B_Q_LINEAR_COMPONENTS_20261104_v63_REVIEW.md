# Qの五つの線形閉成分の独立有限代数監査

2026-11-04。対象は `Q_EXCEPTION_FACTOR`、`Q_FIRST_COMPONENT_UNIT`、`Q_LINEAR_BATCH` の凍結packetである。前の有限例外reviewへhashを結ぶ。著者solver、追加点、素数、β probe、GBは実行しない。

## 完全なsource因子積

factor入力が保護された同じsafe SF37であることを比較する。全九因子のQQ content・指数・全係数積を独立に再構成し、元37次式へ戻ることを確認する。次数目録は1,1,1,1,1,6,6,9,11である。SF37のsquarefree性は既受理なので、この積に重複や未保存の成分を落とす余地はない。

最初の五線形因子から、閉idealの値e=−1,2/3,−1/6,−5/9,−8/9を直接導く。任意の数値点を選んだのではない。残り四因子の既約性を新たに探索する必要はなく、その全積32次を未処理coverとして保持する。

## 各閉fiber全体のsource計算

各線形idealを元のQ(a,e)へ適用し、aモニック二次式を全係数で構成する。e=−1ではQ=a(a+26/225)となる。各保存clopen除去について、除去因子×残り=元modulus、および除去因子がaを割ることを確認する。残りmodulusとaのglobal Bézout=1と、その有限代数内のa逆元も確認する。従って除くのはFULLU外のa=0部分だけである。

残るaの次数は1,1,1,1,2である。最後の二次a代数も両根とnilpotentを含む商環のまま扱い、一根の標本へ置き換えない。各base A_a上へ元a^5Hを直接写し、a^-5で正規化してdモニック四次を構成する。全係数と保存Hを比較する。基底a^i d^jで次元Nは順に4,4,4,4,8となる。

独立QQ quotient-algebra実装で元の全22minorとUを評価する。eの冪は閉idealから導いた有理数の冪であり、generic field root式を例外点へ代入する操作ではない。全22N ideal列をsource image×全基底から再構成し、保存行列の全係数を照合する。さらにUをN回掛けた値と、元22sourceに対する保存cofactorのdot productを確認し、各成分で

`sum_i C_i P_i = U^N`

を得る。membership solverや新しいpivotを探していない。Uを反転すると単位idealとなるため、五つの閉fiber全体のFULLU部分にactual B-bad点は無い。これはgeneric一点を確認した結果ではない。

## 保持する残余

L全成分は前reviewで除去済みである。今回の五線形成分を除いて、Q側の全残余積はe次数6+6+9+11=32。Qがaモニック2次、a≠0でHがdモニック4次なので幾何学的候補は高々256点となる。旧有限残余の高々1600点を保持し、全FULLU B-bad点数上界は**1856**となる。

高次四因子のfiberと旧残余は未判定である。空性、scheme length、全域B rank3、p=0のraw H7、full cubic57、quartic239、Poisson解尽くしへの昇格はない。全input/evidence hash、各閉fiberの検査結果と実測値は同名REVIEWへ保存する。
