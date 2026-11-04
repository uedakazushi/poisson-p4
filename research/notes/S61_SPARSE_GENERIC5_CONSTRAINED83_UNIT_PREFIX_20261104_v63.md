# constrained83 行列の定数単元 prefix と厳密な残余

2026-11-04。入力は凍結済み実236×83行列、49039個の QQ 五変数係数項である。
source C,d,N,L の元配列や T/P、旧692953項 selector を native worker に保持・代入
しない。native が読むのはこの実行列だけで、密236×236逆行列も作らない。

## 計算前に固定した選択と上限

各段で literal QQ 定数 entry があればそれを優先し、無い場合だけ parameter
Laurent 単項式を許す。候補順は、消去で新たに埋まる行列位置の数、項積の上界、
係数 height、元 source 行番号、元 source 列番号による辞書順で固定した。
単項式 support の matching を pivot の存在証明に使ってはいない。

全 pivot の実行前に、その全行更新を support 集合で予測する。product の
Minkowski support が k 項、old entry との和集合が b 項、現在から予測した
所有係数数が L なら、L+2k+3b+16 を scratch 上界とする。これは product、
Laurent shift を揃える二つの作業多項式、結果、および back-add 恒等式確認を
含む。support 集合自体の作業領域も150000項に照合する。

係数の有理数高さ・符号付き exponent の文字数から、次の全 receipt と全残余
を直列化した byte 上界も計算する。既存出力と1MiB metadata reserve を含め
32MiB以内であることを要求した。whole pivot の事前判定が落ちた場合は、
その pivot の一部分も実行しない。native47秒で新しい pivot を止め、最終保存を
含め hard60秒・2GBとした。計算前の policy/input は preflight.json に固定した。

## 完了した実 prefix

一回の native worker は0.414秒で QQ 定数 pivot を24回実行した。
parameter 単項式 pivot は0回で、parameter の新たな分母は生じていない。
全23742行更新について new+factor*normalizedPivot=old を QQ 多項式として
直接照合した。

元行順は

```
188,154,134,157,149,139,159,172,136,185,208,200,
190,210,223,187,202,151,161,0,212,19,32,12
```

対応する元列順は

```
37,17,3,18,8,1,11,9,14,23,38,28,
21,31,29,34,32,12,13,63,33,82,52,53
```

である。pivot の正規化行と全 target factor、各更新の support/height/byte
forecast は一段ずつ receipt として保存した。alpha/beta に課した L 二組も
同じ行操作を受けて保持され、任意の制約行を捨てたわけではない。

実際に所有した係数数の最大は125155、完了 pivot の事前最大は125333。
残余は厳密に212×59、6657非零成分、124840項、最大一成分175項である。
56個の零行も削除せず保存した。列数を減らした一方、総項数は49039から
124840へ増えたことを明記する。残余 rank59は既知rank83と24単元対からの
帰結であって、新しい rank/minor 検査ではない。

## 正しく停止した候補

次の決定的候補は元行23・元列65の QQ 定数 pivot だった。係数 scratch の
予測146462は150000未満だったが、全出力 byte 上界37830400が32MiBを超えた。
従って実行前に停止した。実際に保存したファイルは約5.04MBだが、保守的な
予測を事後に弱めて再実行したことはない。定数 pivot の最大可能数が24、または
他の順序でも停止するという主張ではない。単項式 phase へも進んでいない。

## 戻し方と残る仕事

receipt の行操作を任意の236成分 RHS に同じ順で適用し、将来正当化される
方法で残余を解いた後、pivot 順を逆にして

```
x[pivotCol]=normalizedRHS−sum(normalizedPivotRow[j]*x[j])
```

と戻せる。これが有限の factored inverse instruction であり、密 left inverse を
展開したという主張ではない。実 raw134 RHS と下102個の零制約に適用すれば、
残余を解けた場合に同じ42 normal 列へ戻る。現時点では59列残余の逆、残り84
coupled 方向、ideal 所属、飽和、残余の空性・一般 fibre の尽くしは未計算である。

非単項式 field pivot、別候補順、上限を拡大した retry は実行していない。
