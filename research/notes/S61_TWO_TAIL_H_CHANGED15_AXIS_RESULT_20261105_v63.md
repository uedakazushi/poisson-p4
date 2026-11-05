# 列構成を変えたgenuine15行minor：固定商環では零余り

2026-11-05。ROOTの全source/code監査
`S61_TWO_TAIL_H_CHANGED15_ROOT_SOURCE_REVIEW_20261105_v63/report.json`、
SHA `95d94bfb53418cc6094e8d2bdacdb704496cacf627d98f7389830007306cae74` 後、
許可された一回だけの値workerを実行した。
全450 Laurent行GL逆の積、current14移動による三交換、全正規化、85保存axis係数を
実行前に確認した監査であり、値の独立peerではない。

一次報告 `S61_TWO_TAIL_H_CHANGED15_AXIS_ENTRY_20261105_v63/report.json`、
SHA `2b93dce6522f3f08300270b278bbb860f0cb29e7cbfcc31c9a777ed4bed275ee`。
同じFp32003[a]/monic(H_axis)で

```
det[anchor22,H27/L,Btail1,Btail2] = [0,0,0,0],
D12*detPair4                    = [0,0,0,0]
```

となった。genuine source恒等式は `A*D15_new=+detPair4*D12`。
同じzero-divisor商環で必ず結合積を調べ、別々の非零性からの推論は行わなかった。

全16subsetと16入力ベクトルを保存した証拠SHAは
`9810d5345f635568c974a16a1bd8dcd8ca22ca0f5dc3618f0c49161450ac5dfa`。
1974 field products、owned2231、30 bit、約0.0076秒、入力込みRAW2513542 bytes。
guard停止・retry・別minor/H31・axis・prime・新source評価はない。

この固定検査はinconclusiveであり、char0のH除算・actual rank低下・bad Fittingを証明しない。
既存H-squarefree/厳密GCD/H∤D12と、選択二tail十分条件の非空gapはそのまま保つ。
全H/FULLUのfull57/239、補集合消滅、一般成分の尽くしは未解決。
新たなminorを自動実行せず、H自体の小さい式について次のsource inventoryを別に行う。
