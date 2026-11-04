# Actual零次数normalizer tensor：free4/48と正当なproduct-conductor比較

2026-11-04。元D0のsocleをnormalizer base changeで消してはならない。
そのkernelを明示的に分けた上で、同じR0上の小さいactual比較対象を構成する。
二族のtrace basisは同一視しない。

## 1. 非負gradingからのactual D0

共通横断ambientのzero-degree環をP0と書く。canonical126のzero36方向を実際の共通zero15群切断で除くと
P0はregular21変数環。九つは共通physical R0座標、残り12は八つのzero odd rowsと四つのzero graph rowsである。
full Dの二つのstable idealはhomogeneousであり、非負gradingなので(I_A+I_B)_0=(I_A)_0+(I_B)_0。
従って

```
D0=A_A,0 completed_tensor_(P0) A_B,0.
```

これはnormalizerへのtensorがflatだと仮定する主張ではない。

## 2. 実際の両rank2 normalizer

既存pure CA48/CB48 unit系は、**同じ九physical R0函数**上で

```
B_A=R0+R0*tauA, tauA²=nuA,
B_B=R0+R0*tauB, tauB²=nuB
```

というwhole regular normalizerを指定する。source、cut、proper finite fiberは既受理gateを使用する。
全部のphysical outputは同じsourceから a_Aj+b_Aj*tauA、a_Bj+b_Bj*tauB と与えられる。
そのため

```
E4=B_A completed_tensor_(P0) B_B
```

は、basis(1,tauA,tauB,tauA*tauB)のfree4 R0-algebraを、十二output差

```
(a_Aj-a_Bj)+b_Aj*tauA-b_Bj*tauB
```

とその四basis倍で割ったactual algebraである。正確に48列のmodule presentationを保存する。
八odd rowsは18,19,28,29,60,67,72,63、四graph rowsは64,68,73,74。
後四本を落とさない。

zero physical raw40座標の残り四つはdivergence4式で上記36から復元できる。
その四次minorと全reconstruction matrixを保存するので、12個だけ使うambient座標上の根拠も有限である。
二つのtau乗法matrixの可換性と各二次式をliteralに検査する。

## 3. D0 free9からのmapと係数函数の意味

actual C0 free9 basisは(1,z7*tauA,gA*tauA)。最後のgAは既存CA tail-prefix0の57-unit系が定める函数である。
これを元の48-variable表現へfreshな非単元除算で戻さず、その同じ57系をそのまま保持する。
したがって係数を定めるunit未知数はCA48+CB48+prefix57=153個、全式を保存する。
mapの九列は1と(0,z_i,0,0)、(0,gA,0,0)。
全612 Q0と28 Koszulがこのmapで零となることは、actual tensor productのuniversal propertyから従う。
ここで旧r/g cofactorを新たな0/0商としてcompositionしたり、640列のfresh reductionを実行したとは主張しない。
pruned係数dependency closureは通常のunit scalar演算だけを含み、旧r/g演算を一つも含まない。

## 4. 比較kernel/cokernelに対する正しいannihilator

A_i,0=R0+I_i*tauiであるから、任意のf in I_iに対して乗法はP0-linear map B_i→A_i,0を与え、
inclusionとの合成はf倍である。従って任意のtensor base changeの後も、inclusion mapのkernelとcokernelはI_iに消される。
二つのnormalizationを順に行うと

```
I_A I_B * ker(D0→E4)=0,
I_A I_B * coker(D0→E4)=0.
```

第二段のkernelを第一段のkernelへ戻す際に積を使うので、Torを捨てた操作ではない。
この積がT0全体をannihilateするとは主張しない。Rのgeneric点では両conductorはunitで、比較写像は同型である。

normalizationはfinite integral surjectionなので、tensor後もSpec上surjective。
従ってE4とD0のR0-supportは同じV(F) union P_Rである。

## 5. Conductor曲線とFの交わりの小さいgate

各族のzero7を零にした自身の二変数plane上で、cached actual source2と九座標adapterからR0一次mapを読む。
同R0 N56 unit系の一次係数だけを適用し、v_N=rho61-a_Nのlinear form ellを得る。
nu_Nとz1は原点で零なので、F=v_N²-nu_N*z1²の二次initialはell²である。
cached g_A3/g_B3にell=0の方向を代入し、非零ならgcd(ell,g3)=1。
その場合(I_i,F)はm-primaryで、二変数initial regular sequenceの長さは6である。

**両側でこのgateが通った場合のみ**、product-conductor-ann比較kernel/cokernelのsupportは有限な原点となる。
Rはconductor曲線の外であり、残るV(F)との交点が原点だけだからである。
shared branchが出ればその方向を保存し、有限性を推論しない。

E4自体のsocle/depth、比較kernelが非零かどうか、元D0のdepthはこの構成だけで決まらない。
新しいphysical source jet、point、prime、GB、高次陰関数計算は行わない。

Code `code/ambient_B18_zero_normalizer_tensor_model_20261104_v63.py`。
Report/folder `AMBIENT_B18_ZERO_NORMALIZER_TENSOR_MODEL_20261104_v63`。
