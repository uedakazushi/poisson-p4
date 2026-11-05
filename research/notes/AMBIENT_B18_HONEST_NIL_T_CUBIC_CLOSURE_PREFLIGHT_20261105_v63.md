# 既知一次・二次initial relationsの三次closure：source-only入力

2026-11-05。Sは8変数である。35一次関係×36二次monomialと3実二次NF関係×8変数、計1284列を600次元degree3へ写す。
このpacketは入力と有限縮約だけを作り、三次階数・GB・新しいsource/jet係数は計算しない。
二次行列の独立査読が閉じた後に別のrank jobを実行する。

## 二次NFのsource

既知linear×variableのspanは174次元で、保存済みΠ6はその商を六coordinate NFへ同定する。
実rank3 witness columns280/570/586のΠ6像を同じcoordinateへ埋めた三つのNF quadraticも実initial relationsである。
実sourceの二次項との差はΠ6のkernel、つまり既知の一次関係×変数に属するためである。
三つのpivot NF座標0/1/3を除くと残るNFは`target0*r²,target3*r²,target4*r²`。
Π3とそのunit embedding U3についてΠ3U3=I、kernelは既知全二次span177である。

## 600行を24行へ正確に縮約する

degree3各monomialからsorted順の最初の二変数を因子に選び、Π3で置き換え、残りの変数を付けて
`q:F3→V3⊗S1`を定める。右辺は24次元である。自然な乗法embeddingをEとする。
各`e−Eq(e)`は、二次関係を残りの変数で掛けたものなので既知K3に属する。
従って完全列の代わりに次の有限商の同型を直接得る：

`F600/K3 ≅ F24 / (qK3 + (I24−qE)F24)`。

1284元のq像と24のconsistency relationを全て保持するため、単にmonomial代表を選んだだけの近似ではない。
これは既知closure全体の完全なpresentationである。追加の実三次initial relationが無いことは主張しない。

## 実行前のbudgetと次の判定

原入力は7164 QQ termsだけである。最初のsource-count assertionは三二次式の非零項を11と数えたため停止した。
実際は9項であり、rankや新係数の計算前のcountだけを訂正した。失敗receiptも同folderに保持する。
24行の整数row-scaled入力、row denominator、射影と全source labelsを保存する。
将来のrank workerは24行でstream Gaussianを行い、全列factorizationはhashへ流し、dense24×1308の証明書を増設しない。
固定入力の整数bit boundからminor/inverseのRAW上界も測定する。
得られる`24−rank`は実HF3の上界であり、正確値とは限らない。

既証length16とHF553の下で、残るtailの可能性は`3`、`2,1`、`1,1,1`の三つである。
HF3が零ならNakayamaでtail全体が零になる。HF3=1ならunique standard monomialの因子条件により以後も高々1。
従ってknown upperが1なら最後の形が確定するが、まだその階数は計算していない。
