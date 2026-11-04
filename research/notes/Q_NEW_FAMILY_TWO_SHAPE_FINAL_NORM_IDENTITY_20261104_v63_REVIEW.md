# 既存PRS prefixからの別方式・最終norm identity audit

2026-11-04。元の一回のPRSは、C2/C1を全係数保存した後、2GB制限下で最終段階の
追加メモリ確保に失敗した。その失敗記録は変更しない。
明示的な追加許可によるこの別ジョブは、120秒/6GB/32MiBを上限とし、保存済みprimitive
二次式・一次式**だけ**を読む。PRS全体やresultantは再計算しない。

```text
C2=F2*(a2*s²+a1*s+a0), F2=2^12*r^4*(r−1)^2*T^3,
C1=F1*(b1*s+b0),       F1=2^32*r^11*(r−1)^6*T^8,
Norm=a2*b0²−a1*b0*b1+a0*b1²,
C0=(F1²/F2)*Norm/a2²,
F1²/F2=2^52*r^18*(r−1)^10*T^13.
```

三つの積を順番に足し引きし、square/mixedの一時配列を段階ごとに解放する。
Normをa2²で係数多項式として完全除算し、再乗算で全係数一致を確認する。
既知contentを戻したC0について、元の47,114全整数係数が `raw resultant=−C0` と
一致するかを直接検査する。小さいmonic/nonmonic二例の8×8 determinantの符号も再確認する。

全計算はreversed LEXの5変数係数環で行い、途中の項数・次数・bit数・SHA256と
最終quotient/C0を保存する。多項式gcd、別prime、点、source26、resultant呼出し、
自動延長・別方式への再試行は使わない。完全除算失敗・資源上限ならその状態で停止し、
独立identityを受理しない。実行結果はreportに記録する。

このidentityが受理されても、対象はr射影のraw resultantである。
残る境界因子、rの回収、actual target函数体の指数1、NEWの定義式は別問題である。
既知3/5 degreeの独立証明によるs回収と区別する。
