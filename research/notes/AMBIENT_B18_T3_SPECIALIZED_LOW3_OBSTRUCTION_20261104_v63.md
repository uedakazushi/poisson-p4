# Actual CA平面のLOW3負membership試験：top12は三次までも零

2026-11-04。LOW2でtop12が二次まで全て零だったことを受けて、別に許可された
「最初に検出し得る三次」の試験を行った。保存済みcombined105のsource三次jetのみを使い、
新しいIFT/Newton solveも四次jetも行っていない。

## 有効なのは負の向きだけ

actual共通physical R0からCA zero7平面へのspecializationを行う。
specialization後のQ[3]像に入らないtopを示せば、元のR0のQ[3]にも入らない。
そのようなtopがあれば、既証のT3有限長と合わせてT3非零・cone depth19が従う。
逆向きの推論はしない。とくにTorを無視したT3自体のbase changeとは呼ばない。

三次までの二行jetは20次元の有理vectorである。
Q[3]2034列の全てと、その全許容(r,s)-多項式倍を三次まで含める。
未計算のcoordinate商とg商の微分は独立な正orderの形式的偏差として保持する。
この偏差を任意の正order(r,s)-seriesへ代入した像を含むように、不確定なmonomialごとに
vector-spaceを拡大する。この操作は実際の関係像を拡大するだけであり、
拡大像にも入らないことだけが有効な負の証明となる。

## g_CB tailの三次は実際の座標で結合した

元のfull DAG中のg_CB tailはCBのzero7 inverseの制限である。
その三次式は既受理の

```
g_CB(r_B,s_B)=-6(r_B-s_B)(2r_B+3s_B)^2+O((r_B,s_B)^4).
```

ここへ同じphysical CA平面上のactual CB positive coefficients
r_B=j_CB46,s_B=j_CB95の保存済みjetを代入した。
order3なので、それらの一次部分だけが三次に寄与する。
両者の元node、全保存jet、代入後の三次多項式を証明書へ保存した。
CBのr_B,s_BをCAのr,sと無断同一視していない。
g=rA+sBに使ったA,Bにはorder≥2を使い、他の正則商の偏差はorder≥1とした。
g商の既知原点値0を、微分も0であるという仮定へ置き換えていない。

## 結果と停止点

2034列と全三次までの倍から作る拡大spanのrankは11。
しかしtop12の二行は、未知の偏差を独立変数として保った形式多項式の段階で、
全て三次まで零だった。
従ってこの試験では負のmembership証明を一つも得られなかった。
T3の零性、actual membership、cone depth20やnilradicalについての正の結論は無い。

今回の次数増加は、新しいsource jetを作る探索ではなく、LOW2が消えることと既保存の
source三次解があることを理由に許可された別試験である。
三次でも全て消えたため、この方法はここで停止する。四次や高次を自動実行しない。
残る課題はactual有限module presentationまたは別の正確なmembership証明である。

## 保存・実測

`code/ambient_B18_T3_specialized_low3_obstruction_20261104_v63.py`と
`reports/AMBIENT_B18_T3_SPECIALIZED_LOW3_OBSTRUCTION_20261104_v63/`に保存。
証明書は全Q列の既知vector・偏差の拡大vector、span basisと出所、全topの形式多項式と
残差、g_CB tailのactual bindingを含む。
一回1.256015秒、最大RSS148588KiB。最初の構築時に局所loop名がdegree関数を隠す
実装ミスを検出し、span計算前に修正した。計算次数・入力点・数学的方法は変更していない。
GB、新点、追加のsource計算、p逆元は無い。
