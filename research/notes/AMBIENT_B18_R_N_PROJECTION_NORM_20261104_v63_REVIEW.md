# Actual R0上のN projection norm：独立監査

対象は `AMBIENT_B18_R_N_PROJECTION_NORM_20261104_v63`。actual R曲線とG19切断から共通physical R0へ戻す別の独立レビューを必須依存とする。物理的にRがNの外にあることだけから、R0上の射影像がNの射影像を外れるとは推論しない。

## 保存された二つの単位系を二次まで解釈する

独自の疎 `QQ[T]/(T³)` 演算で、保存済みR19 group系とN56 coefficient系を再評価する。単位逆元は有限geometric seriesを使う。著者workerや単位系のsolverは実行せず、保存された補正を各段で代入して検査する。

すべての変数を明示的にbindし、未指定変数を零にするdefaultは許さない。原点での19×19および56×56 Jacobianを元DAGの独立ADで再構築し、保存inverseの左右積も確認する。次数1・2それぞれのraw residual、補正、補正後の全式三係数を元nodeから照合する。新しい物理source jetや次数3以上の係数は作らない。

R0の九rowは `27,108,109,112,113,117,118,119,61`。NのCohen8は、その同じreadbackの `108,109,112,113,117,118,119,27` である。この順序と実際のseries代入を照合する。

## ノルムの正確な最初の非零係数

Nのactual同R0表示は `tau²=nu, rho61=a61+z1*tau, v=rho61-a61` である。二次までの検算は

```
nu=256T²+O(T³), z1=-23T+O(T²),
a61=-472T²+O(T³), rho61=-240T²+O(T³),
v=232T²+O(T³)
```

を与える。nuとvの定数・一次、z1の定数が実際に零なので、未計算の高次係数はノルム `fN=v²-nu*z1²` のT⁴係数へ寄与しない。従って

```
[T⁴]fN=232²-256*23²=-81600 !=0
```

は有限次数の近似からの推測ではなく、全seriesの正確な非零証明である。

## 同じbaseでのsupportと一般長

受理済みN8 moduleは `R8+m8*tau`、同じR0は `R8[[v]]` で `v=z1*tau`。そのmonic quadratic relationはfNである。一般係数体ではnuが平方でなく、このdegree2 relationが最小であるから、NのR0射影のsupportは正確にV(fN)である。あるいはmonic Weierstrass除算と既受理のactual quadratic normalizerへの代入で、同じkernelを得る。

D4の被約supportはactual N tailのみで、追加nil moduleのcurve supportもN内にある。この事実はliteral W1 tail bindingから使い、全ringをnil支持のcurveだけへ置き換えない。従って、上のfN非零はactual R0曲線のgeneric pointがSupp(D4)を外れることを示す。

T1=T2=T3=0とdegree4安定化は独立受理済みなので、そのR0 generic pointで全positive degree pieceが零になる。ここでは先に有限w-polynomial graded algebraを用い、その後にw=0で完成する。w=0完成環をw=1へ代入する手順は用いない。局所化後は全ringがD0となり、D4への射のkernel T0もD0そのものである。

実際のR曲線のR0像は、別のR-bindingレビューで一つのreadbackがDVR parameterとなり、同じbranchを閉埋込みすることが証明される。この同じbaseで、受理済みの全R\N ordinary-double局所型から、D0とT0のR一般長がともに2となる。nil方向のtorus characterや別baseでの数値を仮定して長さを移したものではない。

T0の平坦性、R全近傍での表示、全support、finite socle、元交叉環の深さは本結論に含めない。必須R-binding gate・全依存hash・実測時間・メモリを専用review reportに保存する。
