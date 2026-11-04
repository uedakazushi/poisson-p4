# Actual D0内の有限socle候補：raw normalizer kernelのTor2長さ63像

2026-11-04。これはproof/size preflightであり、像の非零性はまだ検査しない。
曲線支持の比較kernel全体を有限だとは仮定せず、初めから有限な自然submoduleを取る。

## 1. Physical quotientの前のraw tensor

同じactual R0上でA=R0+I_A*tauA、B=R0+I_B*tauBとする。
先にR0上だけでtensorすると、有限R0-moduleとして

```
Traw=A tensor_R0 B
 =R0 direct_sum I_A*tauA direct_sum I_B*tauB
   direct_sum (I_A tensor I_B)*tauA*tauB.
```

二つのnormalizerのR0-tensorへのmapは、最後の成分だけで
I_A tensor I_B→I_A I_Bを取る。他の三成分は単射。
exact sequences 0→I_A→R0→R0/I_A→0、0→I_B→R0→R0/I_B→0から

```
ker(Traw→B_A tensor_R0 B_B)
 =ker(I_A tensor I_B→I_A I_B)
 =Tor1(R0/I_A,I_B)=Tor2(R0/I_A,R0/I_B).
```

Tor indexは2である。前packetのexact Tor式より、これは

```
A3^15 direct_sum omega_A3^6,
length63, socle36, minimal generator数27.
```

I_A+I_B=(z7,(r,s)^2)に消され、従ってm_R0²に消される。

## 2. 実際のD0への自然map

actual zero ambient P0上のtensorは、Trawを十二physical target差で割ったものなので、自然な全射Traw→D0がある。
よって上のfinite kernelの像W⊂D0は、全次数でfiniteかつm_R0²-annihilatedである。
raw kernelがphysical quotient後にも単射だとは仮定しない。

Wが非零と一つでも検証できればH0_m(D0)が非零になる。
既受理のH0(D0)=H0(T0)とG19 productから、actual Dのcone depthは19に決定する。
逆にW=0でも、他のhidden finite部分が無いとは従わない。

## 3. 二十七代表は八つの一次Koszul cycleだけで指定できる

actual generatorをf=(f1,...,f8)、g=(g1,...,g8)、順序を(f,g)とする。
combined16→(z7,r²,rs,s²,0^6)のunit変換で、六つのzero-column cycle q_alphaを取る。
quadratic3部分のcycle eta_u=(-s,r,0)、eta_v=(0,-s,r)も同じ変換で持ち上げる。
Tor2の27 minimal generatorsは

```
q_alpha wedge q_beta (15本), eta_u wedge q_alpha, eta_v wedge q_alpha (12本).
```

QQ basis63は前15本の1,r,s倍と、後者のeta_u、eta_v、s*eta_uの各六本である。
実際のcycle q=(qA,qB)、p=(pA,pB)に対し、f-blockへの射影からdimension shiftを取ると、代表は

```
sum_i,beta (qAi*pBbeta-pAi*qBbeta) * (f_i tensor g_beta).
```

積へ送るとfull syzygy f*qA+g*qB=f*pA+g*pB=0によって零となる。
したがってD0での代表は同じ係数で(f_i*tauA)*(g_beta*tauB)を足すもの。

選んだ八physical rowsは両側とも18,19,28,29,60,67,72,63である。
D0ではg_beta*tauB=(aA_beta-aB_beta)+f_beta*tauAなので、CA traceで

```
even =nuA*sum h_i,beta*f_i*f_beta,
odd C0 vector=sum h_i,beta*(aA_beta-aB_beta)*image(f_i*tauA).
```

これはold r/gを新たな0/0商へcompositionする手順ではない。

## 4. Quadratic obstructionのために必要なのはcycleの原点だけ

既受理の共通even tangentより全aA-aBはorder2以上。全f_iは原点で零、nuAも原点で零。
従ってq wedge q代表はC0 moduleのodd成分でorder2以上、even成分でorder3以上となる。
その正確なmodule二次項はh(0)だけで決まる。
eta cycleは原点で零なので、eta wedge qの十二代表はodd order3以上であり、二次obstructionへは寄与しない。

q_alpha(0)はsource first coefficientsと、cached plane restrictionのquadratic coefficientsだけから正確に復元できる。
十本はCA first7とCB rows18,60,67を選び、canonical generator(z7,r²,rs,s²)へのconstant matrixのdetは625/768。
その10 inverseから六つのzero-column originsを保存・照合する。
full all-order transformを既に展開したとは主張しない。

## 5. All-order transformの有限source input

七つのcoordinate-prefix105 unit系で、各CB zero coefficientをz7とplane restrictionへ分ける。
さらに二つのaxis105系(r,0)、(0,s)を加える。plane函数Bがm²に属することは既受理なので、

```
B(r,s)=r²*[B(r,0)/r²]
       +rs*[(B(r,s)-B(r,0)-B(0,s))/(rs)]
       +s²*[B(0,s)/s²]
```

の三商は全て正則であり、原点値はcached quadratic coefficientそのもの。
自分の分母が零となる場所へ商をfreshに代入し直す操作はしない。
CB rows18,60,67の三次matrixはunitなので、r²,rs,s²をその三restrictionの線形結合として解ける。
CA f63もorder3だから同じ二次Hadamardで扱える。

これで十selected generatorsのunit matrix、六redundant cycle、二eta cycleがactual source-bound functionsとして指定できる。
既存prototypeは59232 DAG nodes。九cloneの保守的上界は533088nodesで、さらにC0 adapter/mapが必要。
32MiB・2GB内に収まる見込みだが、literal構成workerはこのpreflightでは走らせていない。

## 6. まだ必要なnonmembership gate

代表の有限annihilatorがあるだけでは、Wの非零性は分からない。
最初の候補は十五q-wedge-q代表の二次項を、D0の全640関係の二次initial上界と比較すること。
全linear cancellation carriesを保持し、未計算のregular quotient二次微分をopaqueな自由係数として残した保守的上界が必要である。
それで非所属が一つ出れば実際の非零像となるが、所属してもW=0の証明にはならない。
source5、blindな高次jet、GB、new point、Inv計算は行わない。

Code `code/ambient_B18_finite_Tor63_image_preflight_20261104_v63.py`。
Report/folder `AMBIENT_B18_FINITE_TOR63_IMAGE_PREFLIGHT_20261104_v63`。
