# T1のC四乗商：未実行workerと全multiply-back

2026-11-05。独立Horner peerで確定したv_C(T1)=4を前提に、商workerとliteral入力を準備した。
compileとhelpだけを実行し、商係数・R/C power・multiplybackは未実行である。
fresh ROOT/source-code監査と一回GOを実行前に要する。

入力は既存T1/C Hasse custodyの7595整数／一個のDのbodyをbyte同一で複写し、headerだけ変える。
U1=T1/a³、Q=U1/C⁴、R=−18e²+9e+2。
`(100a−R)^(-4)`の有限Laurent polynomial partから、i=0..16の

```
Qint_i=sum_{j=i+4..20} I_j*binom(j-i-1,3)*R^(j-i-4)*100^(20-j+i),
Q=Qint/(D1*100^20)
```

を使う。負のa冪remainderを「無い」と仮定して捨てない。
完成後に**全係数**で`C⁴ Qint=100²⁰ U1int`を確認するため、負冪に由来する残余も全て排除される。
各a slabは途中e≤73,d≤18を保持し、全和が終わる前にe上界33で切らない。
完成slabの想定Q支持外の非零は削除せずSTOPとする。

source制限Qはa≤16,e≤33,d≤18,全次数≤50の10013 slots、multiplybackのremainderは
a≤20,e≤41,d≤18,全次数≤58の14991 slots。
formulaではinput7595、Rcache289、Q10013、slab1406、metadata16000により54606 cells。
slabを空にして全Rcacheを解放してからC powerとremainderを作る。
検算側の保守的cache289を含む全owned上界は81776で、coefとpacked key双方を数える。

全scalar/root/polynomial積上界2712113。中間整数は全て518 bits以内。
`||Qint||1<4*100^16||I||1`よりQ472 bits、
`||100²⁰I||1+||C⁴Qint||1<2^502`よりmultiplybackの部分和502 bitsである。
したがって518はquotientだけの境界ではない。
45秒soft/60秒hard・2GB・RAW10.61MBを固定する。

停止前にcommitせず、input、全cache、Q、slab、remainder、全cursorをT1C4QP001に保存する。
partial outputから商を主張しない。最終出力はwhole multiplyback零の後だけであり、
出力I/O失敗を第二のfull prefixへ重複保存しない。
固定100 powers/binomialはliteral metadataから復元可能で、未知sourceを再計算しない。

この準備はGCD、他の因子・tail、T0、bad locusの除去を含まない。
既知単元a³C⁴を除くだけなのでFULLUでの零集合は変えない。
