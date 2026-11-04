# Full Bianchi 残余4本の cubic pivots

2026-11-04。前packetの11 quartic pivotsの後に残った4関係だけを扱う。
元319係数DAGは変更せず、448個のgenuine-cubic cofactor coefficientsだけを
exact QQ(a,e,d)で展開した。大きいquartic行列や新点・field・GBは使わない。

## 1. 実結果

前packetの15×15行変換Uの第11–14行は、元blockの行 `0,4,13,14` の
単位ベクトルである。従って残余4はこれら4本のactual319行関係そのもの。
元piece indicesはcertificateに保存した。全4×33=132個のquartic coefficientsは
exactに零である。

残るgenuine-cubic112列、4×112=448個の係数を評価すると **119個が非零**。
既知unit以外を割らない行操作で **4 pivots** が得られた。順に

```
189a²/(2F),  −756a/F,  −36450a/(F*B),  −50000a/(F*B).
```

F、B、Cの意味は前packetと同じである。pivotに必要な逆元はa,F,Bだけ、
元graph/row transportの定義開を含めた新4本の開は **aFBC≠0** である。
unknown factor、d、Jをこの新しい行操作で割らない。

4個のpivot columnsはactual319 source indexで `176,177,237,181`。
source definitionsはそれぞれ

```
G18_normalized*v2, G18_normalized*v3,
G86_normalized*v0, G20_normalized*v0.
```

この4列の正規化済みblockはI4である。4×4行変換をliteral rational
coefficientsとして保存しており、これと元4関係の319係数DAGとの積が
新しい4本のactual全319係数を有限に定義する。残る174個のquadratic-source列のcore係数は
展開しない。

## 2. 独立性と正確な方向数

前packetの新11本には、選択quartic11列でI11がある。そのquartic列で
残余4本は零、上の選択cubic4列でI4を持つ。従って15本のblockは

```
          quartic11  cubic4
new11        I11      *
new4          0       I4
```

であり、未知の `*` を展開する必要なくdeterminant1である。
**直接weight4の15関係は全てaFBC≠0上で点別に独立**となる。

旧quadratic22関係はquartic/cubic両方のsource columnsが零である。
その22本の独立性の既証開はFULLUに限るので、結論は

**FULLU = D(a*d*F*B*C*J) 上で37本の独立kernel directions**。

80方向のgateには **43本不足** する。新15を含むsource-kernel37方向と、
周囲239列への全射性は異なる。actualC122のrank64、全点239下界、
fullquartic/exhaustionを結論しない。

今回cubic coefficientsは零ではなかったので、親課題で挙げられた
「cubic全零ならHq係数も六quad定数独立性により零、残りは旧22に入る」
という代替分岐は適用しない。raw reportのその条件付きfieldとnull fieldは
実際のpositive rank4結果を否定する意味ではない。37方向という正の結論は
上の三角blockから従う。

また、扱ったのは15個の **DIRECT weight4 pieces** であり、lower-weight
Bianchi piecesにnormal変数を掛けた全ての関係、他のsource syzygies、
80次元kernel全体を尽くしたとは主張しない。

## 3. 保存と再現

`code/s61_nonsquare_bianchi_residual_four_20261104_v63.py` を一回実行し、
1.596677301秒で完了した。2GB/60秒/200k DAG node上限内である。
既存10710 nodesのまま、core評価cacheは239 nodes。

certificateは、受理済み親packetのhash、Uの残余4行、元4本のpiece番号、
448個のcubic係数rootと全exact値、132個のquartic零root、4pivotの全行操作、
4×4変換と最終4×112行列を保持する。純粋quadratic/Hq coefficientsは
展開していない。新output先で

```
python code/s61_nonsquare_bianchi_residual_four_20261104_v63.py --output <new-directory>
```

により再現できる。今回二度目のpivot/rank jobは行っていない。
