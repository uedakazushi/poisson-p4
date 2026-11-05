# old25・current15・元full57の大域source橋渡し

2026-11-05。新しい行列係数や行列式は計算しない。ここでの大域とは
`R=QQ[a,e,d,1/(adFBCJ)]` 上であり、`T` や選択minorをRの単元に加えた意味ではない。
既存の旧Omegaに関するファイルは変更しない。本証明の独立監査前には、補集合の削除を採用しない。

## 係数環と元の全source

```
F=500a+90e²+3e−24,
B=(50a+18e²−9e−2)d+60e−40,
C=100a+18e²−9e−2,   J=5ad−6e−4.
```

元63行は六二次式の七light倍42行と真正三次21行である。独立受理済み
ROW_BASISの六関係の削除係数は順に `1,1,1,a,-a,12a/5`、三角minorは
`-12a³/5`。全91列の恒等式なので、残る41 productと真正16行はR上で全63行を生成する。
ここで自由row moduleから全点rank下界を推論しない。

**product rank41は別の全点定理である。** PRODUCT41_FULL_OPENは、独立受理済み
PAIR_EMPTY_QQによる全FULLUのactual RAW24全射と、RAW24_PRODUCT41_BRIDGEの点別橋渡しを
literal hashで結ぶ。後者の六式変換detは80a、ell座標の単元は`-16CJ/243`、
ZWを法とする四式のpivotは`5a,aC/240,16C²/81,160aC/81`。
制限kernelの17方向とRAW24の24方向により下界41、係数10aの非零syzygyにより上界41。
従って**各FULLU幾何点で**元42 productのrankは41である。最後のproduct削除係数12a/5により
選択41行も同じspanである。generic rankという名称だけを前提にしない。

## 元91列からold25へ送る線形写像phi

FULL57_SOURCE_LIFTのphysical91列順と元57行へのLift25をそのまま用いる。
まず`Z=3v0-v1+v6=W=v3=0`とし、五lightを
`(L,x0,y1,y2,y3)`へ変える。逆座標の分母は`-40C`。
実sourceは

```
q1=L²/8100,
q3=L(eL-400a*x0)/8100,
L*x0=(8100e*q1-8100q3)/(400a).
```

従って`(L²,L*x0)`の三次九monomial列を除ける。q2とq4は、既受理の80a変換恒等式により
ZWを法としてq1,q3の組合せである。残るproductの像は、q0と1600q5の五light倍である
Q10のspanにすべて入る。ここでは固定19minorを反転しない。もしold25のrankが25なら、
Q10自身のrank10も同時に従う。

残った26 pure列とH五列に、元のHZ/HWから作る一列
`psi=-wdelta*HZ+zdelta*HW`を加え、phiを定義する。
すべてのproductはHを持たないためpsi像は零である。
genuine16のclearingは`mu_i=-64000*C^(3-cancel_i)`であり全てRの単元。
実sourceで`mu_1=mu_6=delta=-64000C²`、

```
Cdelta=delta*(full57 row47-full57 row42)
```

はphiのkernelに入る。これはfull sourceで零という意味ではない。
残るgenuine15とQ10のphi像は、actual old25×32そのものである。

Cdeltaの元のH pairは `(zdelta,wdelta)`。独立受理済みのliteral恒等式

```
cZ*zdelta+cW*wdelta=eta=192*a*F*B²*C³/125
```

により、この像はFULLU上でprimitiveでありproduct spanには入らない。
zdeltaやwdelta各々、g1/g2/h等の未知因子は反転しない。

任意のFULLU幾何点でold25の25行が独立とする。保持genuine15の組合せがproduct spanに入ると、
phiでQ10との25行関係となるので全genuine係数が零。さらにCdeltaも加えた関係へ同じ議論を
適用すると残るCdelta係数もeta-unitにより零となる。従ってproductを法とするgenuine16は独立。
**全点product41と合わせてrank57以上、元の六source関係により57以下**となる。
これは十分条件であり、actual bad Fitting idealとの等号や逆命題ではない。

## old25からcurrent15への十個の大域unit操作

最初の八pivotは、保存順に

```
-5/4, -5/4, 5a, -5/4, 25a/4, -a/6, -a/6, 2FB²C²/15625.
```

source `(row,column)` は
`(8,20),(1,10),(2,16),(3,12),(0,0),(6,23),(7,24),(11,26)`。
各stepはpivot行をuで割り、保存されたpre-elimination係数ciで他行からci倍を引く。
逆は同じci倍を加えてからpivot行をu倍する。除くのは完全なunit blockだけであり、
全32列のtailを無視しない。pivot順に並べた8×8 blockは対角1の上三角形なので、その逆は
strict上三角部分Uに対する有限和`I-U+...+(-U)^7`であり新しい逆元はない。

第九pivotは固定六係数とpair dualの内積1を使う。dualに現れる分母は非零QQとaの高々二乗。
実GL6のdetは-1、逆は保存済み`E^-1 U0`。pure target `[1,2,3]`の三tail操作とその逆を
含め、全17×24に作用させる。変換後unit行・列を削除して16×23を得る。

第十pivotは実三係数b0,b1,b2とdual`((108e+42)/(5a),-18,-9)`の内積1を使う。
GL3のdetは-9であり、逆を明示した上で**十五pure tail**を全て消す。
全16×23に作用させ、unit行・列の削除後のactual15×22を得る。
二つのsource clearingは`I+ue0^t`（u0=0）だから逆は`I-ue0^t`。
五heavy quotient列`[27,28,29,30,31]`は二completionで変えない。
これらはただの部分行列選択ではなく、すべての行を運んだ3914-nodeの有限DAGである。
その演算は和・差・積のみで、未知係数の逆演算は存在しない。

従ってRの各剰余体で厳密に

```
rank(old25)=10+rank(current15).
```

Q9の全22成分、選択九pure列以外の八tail、H3三列、tail2は保持される。
個々のpure kernelを元rawH7に持ち上げる際は、保存済み
`(cZ*HZ+cW*HW)/eta`倍のCdeltaを引く。この操作の純三次tailへの効果も省略しない。
raw63のQQ正規化と原physical91列順もFULL57_SOURCE_LIFTの受理済み入力に固定する。

## direct M12との結合は別の新chartとしてのみ扱う

DIRECT_M12_FULL15_SOURCEの独立peerは、係数環Rで

```
A*D15_normalized=-T*D12
```

を証明した。H2/H3のLとgamma、a、全非零QQ row倍率、pure tail2のBだけを用いる。
`L=aFC², gamma=-4FBC²/625`。Aの逆を仮定しない多項式恒等式であり、
κ3やP_K92の逆も入れない。新しいGLOBAL橋渡しの独立監査後には、
FULLUの`D(T*D12)`上で右辺が単元なのでAとD15も単元、current15 rank15、
old25 rank25、元full57という十分条件が得られる。

この新source packet自体は、その独立監査を完了したとは記さない。
旧Omegaの条件付きファイルを変更せず、Delta407/Saux/N6/N7barの削除は独立監査後の
別corollaryに限る。T=0、D12=0およびFULLU外の補集合、全FULLU57、全quartic239、
generic五点尽くし、全成分上界は残す。新しい行列値・det係数・点・標数・GBは計算しない。

再現入口は`code/s61_nonsquare_prePOWER_global_full57_source_20261105_v63.py`。
実source・既受理peerのhash、八pivotのliteral factorization、全completion葉の分母、
forward/inverse operatorと全tail参照を保存する。実行済み係数再計算の代替を称さない。
