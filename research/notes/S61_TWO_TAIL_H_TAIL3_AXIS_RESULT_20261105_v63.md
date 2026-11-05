# 固定第三tailの結合検査：完了した零余り

2026-11-05。ROOTのsource/code監査 `99f11d4863e349fc677ac4c037f0001344414b7eeb338603d0db31b7e92efc84`
後、original32 target3の同じ149-node sourceだけを一回実行した。
報告 `S61_TWO_TAIL_H_TAIL3_AXIS_ENTRY_20261105_v63/report.json`、
SHA `4a35df1237b235b7acf62bc9eeb7595aeb79de171ecf157bc11763792cddc23b`。

同じ既存四次商環S=Fp32003[a]/monic(H_axis)で

```
T3 mod H_axis       =[0,0,0,0],
(D12*T3) mod H_axis =[0,0,0,0]
```

だった。零因子環で別々の非零性を使う誤りを避け、指定通り結合積を実際に計算した。
162 typed node、14752 literal source係数、cachedκ4、B一乗clearing、全tail効果を保持した。
4653 field products、owned6733、30 bit、約0.026秒、入力込みRAW453896 bytes。
guard失敗・retry・別prime/axis/tailはない。
これは固定商環でのinconclusive結果で、H|T3のchar0証明でも、rank低下でもない。

## 既計算stateの構造診断

同じ保存stateのnode151..154にあるcachedκ四係数は、全て `[0,0,0,0]`。
一方、四列source node66,94,122,150はそれぞれ非零の四次商環ベクトルだった。
従ってこの固定商環で同じκを使う他のtail収縮も自動的に零になる。
この診断は今回の一次stateの読取りであり、κがchar0でHを因子に持つとは結論しない。
同じ方式で残りtailを順に試す自動fallbackを行わない。

source-onlyの次候補は、新tail値を増やすことではなく、既計算tail1/tail2を再使用して
pair4の四列を `[anchor22,H27,B*tail1,B*tail2]` に替えることである。
旧 `[H27,H31,anchor22]` 三列のκ消滅という特定の原因を避ける方法である。
旧九pure列のrawpair rank1とH3零を保持し、同じM12へ三追加列を付けるgenuine15minorの
多項式Laplace/sign/known-unit証明を先に固定する必要がある。
結合積D12*detPair4を調べなければ同一H成分のopenは証明できない。
この新minorの値計算は一切実行せず、source/costとfresh一回GOを別に要求する。

V1 preparationではunused literal-oneの一項までcopyしたためsource-count assertionに止まった。
当該native/literal部分、失敗code/note snapshot、failure.jsonは元custody folderに保持した。
V2は正しい149-column closure+cached4だけを取り、unused-oneは非activeとした。
値workerの数学的実行はV2について今回一回だけである。
実bad Fitting、全H/FULLU full57/239、補集合消滅の主張はない。
