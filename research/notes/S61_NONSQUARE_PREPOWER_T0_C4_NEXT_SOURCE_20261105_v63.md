# T0のC四乗後：別five-jet Hornerと整数商のSOURCE

2026-11-05。primary k0..3零／k4非零646項の後のSOURCEのみ。
独立Horner値、商、他因子、GCDは計算していない。

## 独立peerの別算術

original root239.polyのQQ7208係数をD197へ直し、既にa⁴を除いたbodyをa次数へ分類する。
producerのR powers/binomialを使わず、truncated Hornerの五jetを伝播する。
R=2+9e−18e²として、J0=I20、J1..4=0から、i=19..0で

```
r=4,...,1: Jnew_r=R*Jr+100*J(r-1)
Jnew_0=R*J0+100^(20-i)*I_i
```

を用いる。r降順の更新でlowerjetを旧値のまま保持する。
各stageの全jet共通denominatorはD*100^(20-i)、最後はD*100²⁰。
crossjetの100倍率を落とさない。入力7208、五jetと一temp各1539、coef/key双方とmetadata16000で
owned48884。606565 scalar products以下、hard518、RAW4MBの監査sourceである。
J0..3零とJ4の全646係数・denominatorを確認する。高いjetの零は主張しない。

## 次のC四乗商は最初から元のDで整数化できる

Cはprimitive整数多項式であり、独立C-valuation確認後にはGaussから
J0=I0/C⁴が整数多項式になる。有限Laurent polynomial partで一slabずつ

```
S_i=sum_{j>=i+4} I_j*binom(j-i-1,3)*R^(j-i-4)*100^(20-j+i)
J0_i=S_i/100^20
```

を作ればよい。partial slabのe≤72を全て保持してからfinal支持外の零と100²⁰の厳密整除を検査する。
最後に**全係数でC⁴ J0=I0**を確認し、負のa冪remainderも消す。
scalar content探索や未知の逆元は要らない。

最終a16/e32/d18/全次数49で9690 slots、slab1387、remainder14592。
全R cacheを解放してから検算cacheを作ればowned79558、全整数hard518、
all scaled/root/product上界2581132、RAW11MBの新しい予算提案となる。
scaledslab435 bits、整数J0は303 bits、multiplyback partialは333 bits以下である。

元のD0は既存J1のD1と文字通り等しい197-bit整数なので、将来の二tail比較にQQ cross-productを要さない。
T0=a⁴C⁴J0/D、T1=a³C⁴J1/DというknownFULLU単元だけの圧縮になる。
J0はまだ未計算。独立C peerとfresh source/code監査・親GOを実除算の前に要する。
bad locusや補集合の除去、strict cover enlargementは含めない。
