# 一個の15次cofactorによるdet認定と旧残余lengthの独立監査

2026-11-04。元norm16の完全QQ必要cover監査に続き、一個の固定cofactor00だけを独自Bareissで計算する。full16 determinant、author worker、solve、新しい点・素数・GBは再実行しない。

## det認定の追加gate

NはQQ[a]上の16×16行列で、前監査が全係数でN*V=delta*e0を確認済み。今回は行0・列0を除いた15×15 minorを、FLINTのinteger polynomial基本演算だけを使う自作Bareissで求める。最後の非零行をpivotに選び、全exact除算を再乗算して確認する。完成多項式と全source NおよびV0を保存する。

cofactor00=V0!=0が成り立つと、普遍adjugate恒等式から

```
det(N)*V0=delta*cofactor00=delta*V0
```

を得る。QQ[a]はdomainなのでdet(N)=delta。これはNの可逆性や事前のdet認定を前提にしない。

## Smithと二つのcover

全16 monic除算の独立監査により、MはQQ[a][e]/(h1)におけるk1乗法であり、coker(M)=QQ[a,e]/(h1,k1)。h1はe-monic16、Mの全係数には負のa冪がない。N=D*MのDは非零QQ定数だからcokerは同じである。

PID QQ[a]のSmith標準形により、det次数148はこのcokernelのQQ長148を与える。a-primary summandの長さはdetのa-order31。aを逆にした明示projection商のQQ長は正確に117となる。

保存されたHはaが単元の係数環上でd-monic4であり、その唯一のd4項が1であることを実係数表で確認する。したがってこのprojection商にHを加えた**明示old containing cover**は正確に長さ468。Hの負冪は既知単元aだけである。

残りのQ例外はe多項式の次数和32、各e値でa-monic2とd-monic4のcoverを持つので幾何学的点数は高々256。受理済みの他の包含・L/Q linear単元除去と合わせたBbadの幾何学的点数上界は468+256=724となる。724は全Bbad scheme lengthではない。117と468のexact lengthも明示した二つの商代数だけについて述べる。未知点の実在・空性・被約性は結論しない。

一回60秒・2GB・新artifact32MiBの独立cofactor計算と、有限なsource monicity確認だけを行う。

## 確定結果と保存

一回の独立cofactor計算はnative3.3514秒、全体4.3010秒で合格した。全1015回のexact除算を再乗算し、完成cofactor00と元V0の全係数が一致した。従って上記domain議論からdet(N)=deltaを独立に受理し、148・117・468のexact lengthと幾何学的上界724を受理する。

最初にscalar fmpz宣言のheaderを欠いたcompile失敗があり、native実行は始まっていない。そのsource・入力・logを別failure packetに保持し、fmpz.hだけを補って同一SHAの整数入力で上記一回の計算を行った。成功証拠と事前compile失敗は混同しない。
