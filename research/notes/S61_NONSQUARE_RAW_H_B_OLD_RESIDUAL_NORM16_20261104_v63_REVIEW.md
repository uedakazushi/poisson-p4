# 旧有限残余norm16の入力・随伴列・元Bézoutの独立監査

2026-11-04。対象primaryはOLD_RESIDUAL_NORM16_INPUTおよびOLD_RESIDUAL_NORM16_CANONICAL。行列式やsolveを再実行せず、保存された全QQ係数恒等式を別の整数convolution実装で照合する。

## 検査する数学的内容

H=h1/LC(h1)、K=k1/LC(k1)はQQ[a][e]のe-monic16次・25次式である。保存Mの全256成分と全16商q_jについて、

```
K*e^j=H*q_j+sum_i M_ij*e^i, j=0,...,15
```

を確認する。反転するLCは非零QQ定数だけであり、未知a多項式を反転しない。

numeric共通分母DによりN=D*Mを整数多項式行列とし、保存16列Vと非零deltaについて全16成分のN*V=delta*e0を確認する。さらにV(e)=sum V_j e^j、Qsum=sum V_j*q_j、二つのmonic cofactor、元のQQ leading unitを戻したcofactorまで確認し、最終的に

```
delta=original_h1_cofactor*h1+original_k1_cofactor*k1
```

を元sourceの全QQ係数で直接検査する。これだけでdeltaが元idealに属することが分かり、deltaをdeterminantだと認定し直す必要はない。primaryのdet objectとdeltaの保存係数一致は確認するが、det計算の独立再実行と呼ばない。

deltaは非零a単項変数式、次数148、最小a指数31。整数coefficient gcdでQQ contentを確認し、それとa31だけを除くと、constant term非零のprimitive117次多項式が残る。a非零の旧h1=k1=0残余に対する必要coverである。未知a因子を逆にしたり除去したりせず、因数分解もしない。旧残余、他のBbad成分の実在・空性、FULLU全rankは結論しない。

## 独立演算と保存事故の扱い

専用Pythonはauthor helper、native worker、determinant、solveをimport/実行しない。各多項式のQQ分母を整数定数で払い、Python任意精度整数の係数convolutionを行う。a,eを6bit幅で符号化するが、全積でe次数が64未満であることを先に確認するため、指数carryによる同一視は起きない。authorのfmpq_mpoly canonicalityに依存しない。

最初のsort-without-combine失敗はfailure receiptとsourceを保持し、成功証拠として使わない。後の全数学完了後のraw aggregateは34721772 bytesで32MiB guardに失敗しており、その失敗も保持する。lossless gzip公開は計算の再実行ではない。両gzipをstream復元し、raw byte数・SHA、保存済みlocal rawのSHAを照合する。これにより保存修正と数学的再計算を区別する。

一回60秒・2GB・新artifact32MiB上限。新しいβ点、素数、GB、16次detの再計算は使わない。

## 確定した独立検査結果

一回4.1424秒、最大RSS79320KiBで全検査に合格した。独自整数係数積1961810回、検査中の最大support5094項。全16除算・16個のNV成分・元Bézout・QQ content/a31・lossless復元を受理する。ここまでの独立結論は非零117次必要coverであり、determinantの独立認定を必要とするexact lengthは別gateに保つ。
