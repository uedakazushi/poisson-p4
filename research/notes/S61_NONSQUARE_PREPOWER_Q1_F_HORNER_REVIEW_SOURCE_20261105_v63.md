# F valuation零の別Horner独立監査SOURCE

2026-11-05。producerのbinomial／R powers engineと異なる算術の監査sourceのみ。
このpacketはHorner値を計算せず、ROOTまたは別peerの実行を要する。

ROOT正規化済み4811整数係数をa次数ごとの17個の(e,d) mapへstream分類する。
R=24−3e−90e²、共通denominator D0*500¹⁶の下で

```
H=I16;
for i=15,...,0: H=R*H+500^(16-i)*I_i.
```

とすれば、Hはproducerのk=0 polynomialそのものである。
原binomial、Hasse derivative、全R powers、native source関数をimportしない。
`||R||1=117<500`から全signed partialはinput_l1*500¹⁶未満、445 bits以内。
input4811と二つのHH map1254ずつ、coef/key双方、metadata16000でowned30638。
全scalar積65020以下、RAW4MB、45秒soft/60秒hard・2GBの別監査予算である。

計算後input mapを解放し、保存済み665係数をstream比較して全てとdenominatorを照合する。
zero-ordersはk=0が最初なので空ファイルだけを確認する。
F因子が無かったため、F除算のsourceやworkerは準備しない。
他因子、GCD、他tail、不良集合への新しい推論も行わない。
