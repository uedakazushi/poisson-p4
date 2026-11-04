# Zero7を保持したhonest9全九変数一次試験：top12は全て一次まで零

2026-11-04。CA平面へのspecializationを止め、同じphysical R0の全九方向を保持した
正確な9行presentationの最初のnegative testを実行した。
zero7方向とKoszul kernelを残すため、以前の二変数試験とは異なる。
T3非零・零、全D depth19/20はいずれも今回決定していない。

## 二つだけのpure57系で全coordinate商の一次係数が決まる

actual CA/CB各自のF9を(Y0,...,Y8)とし、それぞれのsource48式とF9式の計57式を
QQ[[Y9]]で二次まで解いた。各原点Jacobian逆行列は一度だけ使う。
保存したmonomialは定数1、九一次、45二次の計55個。
両方の57式について全55係数を最後に照合した。
八つ、十六個のprefix inverseを別々に展開していない。

pure source functionが

```
b(Y)=sum L_i*Yi + sum_(i<=j) Q_ij*Yi*Yj + O(m³)
```

なら、prefixを一つずつ増すHadamard係数h_iの一次jetは

```
h_i=L_i + sum_(j<i) Q_ji*Yj + Q_ii*Yi + O(m²).
```

r,sを常に残すzero7 telescopeでは、さらにQ_i7*Y7+Q_i8*Y8を加える。
この有限公式から新CA商792個と旧CB coordinate商784個を全て求めた。
各定数項が凍結済みの原点値と一致することも確認した。
g=rA+sBのA,Bはgのorder3から一次まで零である。

CBの一次座標をCAの同じphysical baseへ戻す線形adapterも保存した。
そのdeterminantは **1/4096**。
共通physical x9の一次jetが二経路で一致することを確認したので、
別々のabstract normal formの座標を同一視したのではない。

## 四つのg商の未知微分は実際に消える

四つの旧g商の原点値0は既受理だが、その微分を0とは置かない。
四つの独立なorder1形式偏差として保持した。
元DAGのg tailはorder3なので一次には寄与しない。
全2069 honest9 relationと全top12を、九つのCA座標と四つの偏差の一次式として評価した。

その結果、**四偏差の係数は全relation・全topで全て正確に零**だった。
従って保存した一次matrixは、単なる不確定微分を加えた上界ではなく、
実際のhonest presentationの正確な全一次matrixとなる。
ソース二次より高いjetや、g商を割るための四次jetは使わない。

## 結果

90-dimensional vectorは九generator各々の定数と全九一次係数である。
全2069 relationは原点で零なので、一次像には定数倍だけを含めれば十分。
そのspanのrankは **43**。
しかしtop12は、零weight七方向も含む全九方向に対して、全成分が一次まで零だった。
従ってnegative certificateとなるcolumnは無い。

この結果はT3=0またはmembership成立を意味しない。
T3自体のpresentationは依然必要であり、元free module中の消失orderからT3の長さや
nil指数を推定することもしない。
先行結果のT3有限長<=36は維持するが、その零・非零は未決定。
自動的なhonest二次jetやsource三次への拡張は行わず、この一次方法をここで停止した。

## 保存・計算量

Source `code/ambient_B18_T3_full9_firstjet_20261104_v63.py`。
`reports/AMBIENT_B18_T3_FULL9_FIRSTJET_20261104_v63/full9-firstjet-certificate.json`に
二つの57-system/Jacobian/逆行列/全source2jets、全関数2jets、792+784商の一次係数、
adapter、全2069×90 matrix、四偏差の全感度、span basis、全top12を保存した。
一回7.050835秒、最大RSS500728KiB。
60秒/2GB/32MB内。GB、新点、source二次を超える、またはhonest一次を超える計算、
nonunit inverse、CA平面へのspecializationは無い。
