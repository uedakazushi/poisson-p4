# Co-Higgs rank M三・tau零：全二源を二normalの有限問題へ

2026-11-08、PRIMARY E、log_manuscript。有限job前の全source入口。
1054主稿／frontierとglobalは凍結したままである。

## 1. 元の全源と、m二の古い除外を全mへ戻す根拠

全multiplicityでconormal可換とする既受理local modular補題から、
O(1) co-Higgs sourceの全35/30/9と全Lax式が有効である。
rank M三ではM:C3→sl2は可逆、R=tau M(U)であり、
全actual normalizerによりsemisimple／nilpotentの二whole symbolへ尽くされる。

base(a,b,c)、normal(z,w)、M=[[b,c],[a,−b]]として、
全weight零sourceは

~~~
P0=T∧(z∂z+w∂w)+U∧((bz+cw)∂z+(az−bw)∂w)+K∂z∧∂w.
~~~

semisimpleはU=∂b、T=2tau(a∂a−c∂c)、K=(1−2tau)zw。
nilpotentはU=∂a、T=tau(c∂b−2b∂a)、K=(tau−1/2)z²。
actual leading五次は−4ac tau(az²−2bzw−cw²)、c² tau(az²−2bzw−cw²)。
従って高multiplicity m≥3は正確にtau=0であり、両型を全部保持する。
旧m二sourceで除外したtau零の保存証明書を、
[all-m独立review](PLANE_HIGH_M_RANK_AT_LEAST_TWO_SEPARATE_REVIEW.md) §7が既に全源で復活している。

## 2. 既存完全kernelから先験的に有限な全模型

保存report CLASSIFICATION_PLANE_RANK_ZERO_FULLRANK_WEIGHTS_20261107.jsonには、
両tau零の完全D1核八、actual shear rank六、D2 rank九が保存され、
上記別作者all-m reviewで元全350/30/9を再構成済みである。

[finite shear補題](PLANE_COHIGGS_SHEAR_EXACT_SOURCE_REDUCTION.md)の同じ実作用を使う。
旧Sの同じ六minorからK=S A+N1 b1+N2 b2を全30×8で読む。
残二列を既存核から取るだけで、新M1、point、rank又はminorを作らない。
全有限shear逆を戻し、任意P1はuN1+vN2となる。
旧D2の同じ九inverseから、任意全Poisson solutionのFULL九P2は必ず

~~~
P2=u²Q20+uvQ11+v²Q02
~~~

に等しい。これはP1についてのweight二Jacobiが二次で、
D2が単射だからであり、formal IFTやjetからの推測ではない。
従って全候補
P(u,v)=P0+uN1+vN2+u²Q20+uvQ11+v²Q02
を元全350式へ戻せば、全未知higher条件が次数高々四の有限二変数idealとなる。
全actual Euler-Pfaffianも全(u,v)の多項式として保存する。

## 3. 成功・非閉鎖の判定を先に固定する

一回のwhole二源batchで、元idealの二次／三次／四次式と全Pfを得る。
原式の手で読める係数組合せがu=v=0、二本のcoordinate axis、
又はproper因子を強制すれば、その全supportを証明する。
all350が全模型上零でも、proper因子又は既知親への実全源mapが無ければ
新しいE成功とはせず、明示した全二parameter sourceを残す。
新point、追加minor、GB又はjet延長へ進まない。
旧D2核が単射と食い違えば全free P2を残して停止する。

予算はwhole二源一回180秒／3GiB。旧caller／旧report／旧jobは変更せず、
既存の標準Fraction sparse primitivesを読む独立inline callerを新job内へ固定する。
ここでの全源候補復元と、残二変数idealの全support証明は区別する。
