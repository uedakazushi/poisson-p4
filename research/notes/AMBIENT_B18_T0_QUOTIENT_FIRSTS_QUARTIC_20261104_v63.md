# D0で必要な商の一次微分：二族の限定quartic contraction

2026-11-04。read-only sensitivity preflightで、D0の640関係のうち四列だけに
`old G_CB,i first - new H_CA,i first`が残ることが分かった。
CA側だけではこの差を決められないため、二族それぞれの保存済みsource3から
必要な四次五monomialだけを一回追加する。
T0のsocle、cyclicity、full D depthをfirstjetから推論しない。

## 実際のcacheと新しい計算量

CA57のsource3は既存joint105 certificateの対応57変数。
CB57のsource3は、旧CB g商原点値packetのsingle actual57 certificateに実際に保存されている。
CBのunknown順序、48方程式、F9順序がFULL9 source2 frameと文字どおり一致することを
quartic計算前に確認した。新しいCB cubic solveやCB48 refitはしない。

旧pure source DAGの必要ancestorはCA9735、CB9978であり、各変数はその族の57個以内。
二変数degree≤4の15monomialを使うが、未知の追加係数は各族57×5=285だけ。
保存source3を四次0で拡張して方程式のquartic残差を測り、保存済み固定J57逆行列で一度補正する。
補正前の全三次方程式と、補正後の全57式・全15係数を確認する。
新しいsource degree5以上やwhole honest二次matrixは作らない。

## 正則商firstのexactな決定

各族でzero7座標を零、positive二座標をr,sとする。
gはf63のこの曲線制限、Bbarはphysical rows64,68,73,74のodd係数の同じ制限。
既証の正則商とH(0)=0から

```
Bbar=g*H,    Bbar_4=g_3*H_1
```

が全係数で従う。g4の寄与はH0=0により無い。
g3のr³/s³係数は非零なので、Bbar4のr⁴/s⁴からH1二係数を求め、
中央三係数も含む全五quartic係数で恒等式を照合する。
これを四つの商、両族で行う。
旧regular-function DAGや原点値は変更せず、exactな一次微分だけを別override証明書に保存する。

CBの二つのparameterをCAの二つと同一視しない。
保存されたactual F9 adapterの最後二行を使い、G_CBのfirstをFULL9 CA座標へ戻す。
H_CAはCA r,s方向のみに依存する。
従って求まる四つの差は九変数全体に対する正確なlinear formsである。

## D0一次関係への代入と範囲

前preflightの既知640×81 matrixへ、四つのold CB商と四つのnew CA商の
正確なfirst derivativeを代入する。全未知感度はあらかじめ保存したので、
source4を使ってwhole scalar二次matrixを再構築する必要はない。

actual firstspan rankと、visible g_A*tau basis classの各一次parameter倍が
そのspanへ入るかを保存する。
この初期moduleの判定だけから、actual elementのannihilator、socle、T0 cyclicity、depthは結論しない。
hidden higher-order torsionは依然別問題である。

Code `code/ambient_B18_T0_quotient_firsts_quartic_20261104_v63.py`。
Report/folder `AMBIENT_B18_T0_QUOTIENT_FIRSTS_QUARTIC_20261104_v63`。
raw source3、quartic残差・補正、全source4、g/Bbar15係数、全商first、全640関係を保存する。
60秒・2GB・32MiBの一回job。点、prime、GB、degree5への自動増加は無い。
