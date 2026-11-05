# T1だけの既知C単元valuation：未実行の固定worker

2026-11-05。T1の7595項を保存済み整数／一個のD表示からbyte同一で切り出した。
四source列を再評価せず、旧149/572sourceも再実行しない。
今回は**source、custody、compile、helpのみ**。C valuationは未実行であり、
fresh ROOT/source-code監査と一回GOを実行前に要する。Tや残る七tailは含めない。

U1=T1/a³はa次数20、e次数41、d次数18。
`C=100a+18e²−9e−2=100(a−a0)`、`a0=(-18e²+9e+2)/100`とする。
U1(a0+z,e,d)の係数をk=0から順に計算し、最初の非零で止める。
Cはaに関しlinearで100はQQ単元なので、そのkが厳密なv_C。
Cはaを割らないからv_C(U1)=v_C(T1)でもある。

common denominatorはD1*100²⁰で、整数演算だけを使う。
各Hasse係数は

```
sum I_ijl * binom(i,k) *100^(20-i+k)
    *(-18e²+9e+2)^(i-k) *e^j*d^l.
```

R(e)=2+9e−18e²のl1は29<100。
入力l1が2³⁶³未満、binomialが2²⁰以下なので、各部分和を含む整数上界は
input_l1*2²⁰*100²⁰、518 bits以内である。
root powerは21段・計441係数、二変数accumulatorはe≤81,d≤18の1558個以下。
入力7595係数とpacked key、root係数と保守的index、accumulator係数とkey、
全metadata/scratch16000を合わせて35188 cellsを維持する。

全21層のforward上界3349395積。最初の非零層だけ逆順に引いて零へ戻し、再現して保存する。
検算二passを足しても3972185積で、先行の653万上界以内である。
45秒soft・60秒hard・2GB、全RAW19.44MBを維持し、guardはcommit前に働く。
停止時は全input、全root slots、accumulatorと全cursorをT1CHP001へ保存する。
部分層やprefixからvaluationを結論せず、自動継続もしない。

ここで商、GCD、他の既知因子F/B/J、他tail、点、標数は扱わない。
最初の非零Hasse層のsource独立監査は、将来の実結果を得てから別途行う。
既知単元valuationだけでは不良集合を除去しない。
