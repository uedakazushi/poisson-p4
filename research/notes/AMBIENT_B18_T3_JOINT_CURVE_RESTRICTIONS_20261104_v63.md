# 同じphysical base上の両zero-odd ideal：最初のcombined105計算

2026-11-04。CAのcurve inverse57系へCBのsource48系を接続し、最初のactual同時restrictionを
三次まで計算した。別々のanalytic normal formを同一視していない。
全Dのdepth、T3の非零性、ell-colonはまだ計算していない。

## Literalな105-unit入力

元のsource DAGをそのまま使う。未知数はCA coefficient48、共通physical zero base9、
CB coefficient48の計105個であり、方程式は

```
CA source48 = 0,
CA (f18,f19,f28,f29,f60,f67,f72,j46,j95)=(0,...,0,r,s),
CB source48 = 0.
```

CBにも同じ九つのphysical baseを渡す。CA/CB各自のtraceを保持し、pは逆にしない。
原点Jacobianはblock triangularで

```
det J_CA57 = -3623878656/625,
det J_CB48 = -331776/625,
det J_105 = det J_CA57 * det J_CB48 != 0.
```

この全105式・全binding・105系の出力としてg_CAとCB zero odd8本を有限DAGに保存した。
固定原点Jacobian逆行列で、`QQ[r,s]/(r,s)^4`の次数1、2、3を順に解いた。
最後に105式の全10係数を照合し、Jacobian逆行列も積が単位行列になることを確認した。
これはこのcombined系の初回計算であり、旧ungauged R9系、高次jet、点、GBは使わない。

## 得られたactual restriction

CA curveの方程式は

```
g_CA = 3r^3/8 + 2r^2s + 13rs^2/6 + 2s^3/3 + O((r,s)^4)
     = (r+4s)(3r+2s)^2/24 + O((r,s)^4).
```

CB zero odd8本をこのactual CA zero7平面に制限したものをB_iとする。
全て定数・一次部分は零。三つの二次初期形は次のとおり。

| CB physical row | 二次初期形 |
|---:|---|
|18|`-95r^2/72 -115rs/36 -5s^2/3`|
|60|`-1385r^2/288 -875rs/72 -145s^2/24`|
|67|`-65r^2/72 -10rs/9`|

row19とrow63は三次から始まり、その全三次係数も保存した。
row28,29,72は三次まで零だが、恒等的な零とは主張しない。
g_CAの三次初期形と全B_iの検出した初期形の共通projective零点は無い。
したがってこれらのactual idealはm-primaryである。これはその初期形が全平面の
原点だけで消えることによるもので、有限jetで全次数の零を推測したものではない。

さらに上記三つの二次初期形は全二次形式を張る。
そのdeterminantとNakayamaによる正確なidealの同定は、再jetを行わない別の有限係数receipt
`AMBIENT_B18_T3_DOUBLE_ANNIHILATOR_LENGTH3_20261104_v63`に記録する。

## Scope

先行packetの両側annihilator証明が与える係数環は、同じphysical R0での
`R0/(I0_CA+I0_CB)`である。今回のliteral105系はその八つのrestrictionを実際に指定する。
今回の計算だけでもこの環が有限長であることが分かり、従ってT3も有限長となる。
ただしT3=0の可能性を排除しておらず、全Dのdepth19という結論はまだ出さない。
ambient二行matrixの原点零性や閉fiberの長さからmembershipを推論していない。

## 保存と実測

Sourceは`code/ambient_B18_T3_joint_curve_restrictions_20261104_v63.py`。
`reports/AMBIENT_B18_T3_JOINT_CURVE_RESTRICTIONS_20261104_v63/`の
`joint105-algebraic-DAG.json`にliteral algebraic入力、`joint-degree3-certificate.json`に
固定Jacobian・逆行列、105未知数の全三次jet、各次数補正、八restrictionとg_CAの係数を保存。
一回1.086132秒、最大RSS94468KiB。60秒/2GB/32MBの制限内。
