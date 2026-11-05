# 既知C四乗を除く有限operator：SOURCEのみ

2026-11-05。実valuation primary v_C(T1)=4の後のsource予算であり、商は未計算。
独立C-valuation peerと新しい親budgetを得る前に係数workerを実行しない。

U1=T1/a³=ΣI_j(e,d)a^j/D1、C=100a−R(e)、R=−18e²+9e+2とする。
Q=U1/C⁴はa次数16、e次数33、d次数18、全次数50。
これは因子C⁴の各次数を引いた積の次数の恒等式であり、partial sumが最初からこの支持に
入るという意味ではない。

`C^(-4)=100^(-4)a^(-4)Σ_h binom(h+3,3)(R/(100a))^h`の有限polynomial partから

```
D1*100^20*Q_i
 =sum_{j=i+4..20} I_j*binom(j-i-1,3)*R^(j-i-4)*100^(20-j+i)
```

を得る。i=0..16の一slabずつ作れば、partial e上界73・d上界18で1406 slots。
完成したslabのfinal bounds外の非零係数は削除せずSTOPとする。
完成Qは10013 slots以下で保持し、R powersを解放してから
`C⁴*Q_integer=100²⁰*U1_integer`を全三変数係数で確認する。
このremultiplicationは25-term以下の固定known polynomialとの積だけで、GCDもfield inverseもない。

input7595、R-cache289、Q10013、検算remainder14991、二つのcoefficient/keyセルと
全metadata16000により最大owned81776。source polynomial productsは2445280対以下。
全ての整数は518 bits以内であり、商l1上界472 bits、検算部分和502 bitsから従う。
`Σbinom(h+3,3)(29/100)^h <(100/71)^4 <4`を用いる。
全RAW15MB・45秒soft/60秒hard・2GBの提案で、まだworkerを構築・実行していない。

これにより次数は減るが、項数が必ず元7595より減るとは保証しない。
FULLUで反転済みのa³C⁴を除くので幾何学的零集合は変えず、補集合は消えない。
次の未知因子や他tailへ自動的に進む計算でもない。
