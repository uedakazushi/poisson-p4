# 確定Hの一回固定六次Rabin検査：実行前source

2026-11-05。本記録は準備のみ。Rabin値・特殊化値・既約性の結果はまだ計算しない。
入力は独立全係数検証済みの94項primitive整数Hだけであり、betaやphysical source oracleを呼ばない。

aの最高次係数は、既存H表により

```
[a^6]H=17578125*d^3*(80-25d+66e*d).
```

ここでL=80-25d+66edはQQ[e,d]で既約。e=0のa定数係数は
`-1025280-12816d` であり、dで割れず、e=0,d=16/5で `-5331456/5`。
したがってLでも割れない。a係数の共通因子は先頭係数のdまたはLを割るはずなので、
HはQQ[e,d][a]でもprimitiveである。

固定する特殊化はただ一つ、e=0,d=1、prime32003。
先頭a係数は17578125*55で、このprimeで非零。これ以外のprime・特殊化は使わない。
HがQQ[a,e,d]で可約なら、primitive性により非自明な二因子はともに正のa次数を持つ。
特殊化後のa次数6が保たれるため、両因子のa先頭係数も残る。
従って固定六次像がFp32003[a]で既約なら、HはQQ[a,e,d]で既約である。

実行時に単項化したfへRabin条件

```
x^(p^6)=x mod f,
gcd(f,x^(p^3)-x)=gcd(f,x^(p^2)-x)=1
```

だけを調べる。六つのFrobenius剰余と二つの全Euclid quotient/remainder chainを保存する。
可約だった場合はこの固定候補のnegative結果で止まり、再特殊化・別prime・factor familyを追加しない。

全多項式は次数6以下。二つの次数5以下の積mod monic fは高々36+30=66 scalar積。
15 bit指数の二進冪は高々30積なので、六Frobeniusは11880 scalar積以下。
正規化と二EuclidのLC逆・division・inverse-product checkを含めても14000 scalar積未満。
全field積と一回加減の一時整数は30 bit以内、固定入力係数は既存40 bit以内。
完全stateとfailure copy・metadataを含むowned上限8192、RAW上限1MiB、45秒soft/60秒hard/2GB。
特定の30 bit上限を既存H入力の40 bitと混同しない。

実行前custodyはsource/worker/このnoteをhashで固定し、--helpだけを確認する。
Rabin値を生成するmode、合成field fixture、一般因数分解を準備段階で走らせない。
実行には別途ROOTのfresh source/code監査と一回GOが必要である。
既約性が成功しても、H=0上の実bad Fittingやfull57/239を意味しない。
