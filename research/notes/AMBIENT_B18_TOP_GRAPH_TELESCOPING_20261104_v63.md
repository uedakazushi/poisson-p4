# Actual cubic graph六pairとtop T3十二pair：単元座標と有限telescopingで構成

2026-11-04。直前のblocked preflightのbytesを保ち、新しい方法で不足したmembership
coefficientsを実際に構成した。大きいGB、source再生成、ungauged R9の再計算、点、
formal-seriesの切断は用いない。完成したのは**係数関数の有限代数回路**であり、
T3への所属・socle・depth19/20・全nilradicalの判定は行っていない。

## 1. 保存済みactual関数の9×10 Jacobian gate

物理baseを

```
R0=QQ[[x27,x108,x109,x112,x113,x117,x118,x119,x61]]
```

とし、CBのactual零odd八本と、positive oddをwで割った二本をF0,...,F9とする。
前packetの48unit Hensel equationsと41,315-node source DAGへ、原点でexact rational
forward ADを行う。48未知数zの微分は

```
dz/dX=−(partial_z H48)^−1 partial_X H48
```

で求め、全10本のactual dF/dXを保存した。source48 Jacobianの実determinantは
`−331776/625`、既知のblock determinantの絶対値と一致する。これは未補正の
source変数をphysical Xと見なす微分ではない。

原点で10本とも零で、Jacobianのrankは**9**。index順

```
0,1,2,3,4,5,6,8,9
```

を選ぶと9次determinantは **1/96**。対応するactual physical odd rowsは

```
18,19,28,29,60,67,72,46,95.
```

除いたrow63の一次部分は零である。従って選んだ九つをY1,...,Y9とすれば、
`Y=F_selected(X)` はpointed formal/algebraic coordinate changeである。
特にactual CB odd ideal `I_B0+J_B0` はR0の最大ideal全体に等しい。
abstract normal formの座標を同じものと仮定したのではなく、actual sourceに束縛した
微分でこのgateを確認した。

## 2. Inverseは九つの有限unit57系で指定する

source48 equationsをH48(X,z)=0と書く。actual inverseを求める有限系は

```
H48(X,z)=0,
F_selected(X,z)=Y.
```

未知(X9,z48)に対するSchur determinantは、z先・X後の順で

```
(−331776/625)*(1/96)=−3456/625 !=0.
```

従って一意なpointed algebraic-power-series inverseがある。
新しい高次seriesを計算せず、同じfinite DAGを九回instantiateし、i番目ではtargetを

```
Y^(i)=(Y1,...,Yi,0,...,0),   i=1,...,9
```

とする。九つの全57式、variable bindings、unit determinantを保存した。
最後にY自身をactual F_selected(X)へ合成するので、全出力は同じphysical R0上の函数である。

## 3. 六つのmembershipを有限telescopingで実際に解く

CBの六つのactual physical cubic odd係数をb_alpha(X)とする。全て原点で零であることも
保存DAGで確認した。B_alpha^(i)=b_alpha(X(Y^(i)))、B_alpha^(0)=0と置く。
次の九係数を**具体的な演算回路**として定義する。

```
h_alpha,i(Y)=(B_alpha^(i)−B_alpha^(i−1))/Yi.     (1)
```

Yi=0に制限すると二つのtarget Y^(i),Y^(i−1)が一致し、unit57系の一意性により
numeratorは恒等的に零になる。従って(1)は正則power seriesであり、algebraic函数の
coordinate quotientなのでalgebraic power seriesでもある。
これはpや未知因子の任意の逆元、積分記号、有限jetからの推測ではない。

DAGではこの演算を通常のunit inverseと区別した `r` nodeとして記録する。
各nodeはnumerator ID、coordinate ID、restriction-zero証明のinstance番号とorigin値を持つ。
有限polynomial equationsへ戻すなら `Yi*h=numerator` とし、上のregularityと
R0がdomainであることから一意な正則解を選ぶ。単元Jacobianを持たないこの最後の式だけを
無条件な任意branchの方程式として使ってはならない。

全54 regular-coordinate quotientについて

```
sum_i Yi*h_alpha,i = B_alpha^(9)−B_alpha^(0) = b_alpha(X)
```

がexact telescopingで成り立つ。最後の等号は、full target Y=F_selected(X)で
元のX,z自身が同じunit57 inverseを満たす一意性による。
選択しなかった10番目のcoefficientを零として、これで実際の1×10 membership六本のsolutionを得た。
**未知hを残したrecipeではなく、六つ全てのhがfinite algebraic circuitで指定された。**

## 4. Same physical base上のCB graphをCAへ代入する

前packetのactual trace basisは、それぞれ

```
tau_A=p_A−s_A/2, tau_B=p_B−s_B/2,
tau_A²=nu_A(X), tau_B²=nu_B(X).
```

二つのtauを同一視しない。raw physical zero odd出力を
`q_i=A_Ai+f_Ai*tau_A` または `A_Bi+f_Bi*tau_B`、positive出力を
`q_e=w*(A_Ae+j_Ae*tau_A)` またはCBの同形と書く。
これら全係数は、同じphysical X上の前packetのactual DAGで既に指定されている。

六つのmembership solutionを使えば、CBのactual raw cubicを持ち上げる一つのgraphは

```
Phi_Balpha=w³*A_Balpha
 +w³*sum_i h_alpha,i*(q_i−A_Bi)
 +w²*sum_e h_alpha,e*(q_e−w*A_Be).
```

CB imageへ代入するとexactにw³(A_Balpha+b_Balpha*tau_B)へ戻る。
これへCAの同じphysical qを代入してCA graphを引いたnormalized差は、**同じtau_A基底で**

```
a_alpha=A_Balpha−A_Aalpha
 +sum_i h_alpha,i*(A_Ai−A_Bi)
 +sum_e h_alpha,e*(A_Ae−A_Be),
b_alpha=−B_Aalpha+sum_i h_alpha,i*f_Ai+sum_e h_alpha,e*j_Ae.
```

六つの(a_alpha,b_alpha)すべてをDAG output IDsとして保存した。最終的な十二pairは

```
(nu_A*j_Ae*b_alpha, j_Ae*a_alpha),
  e=positive physical rows46,95,
  alpha=cubic physical rows9,39,54,84,99,144.
```

これら全12個も実際のfinite algebraic circuitとして保存した。
graph liftの選択差はCB kernelに入るので、全CB kernelとlower graph relationsを保てば
full Q idealは変わらない。それらを削除して十二pairだけをT3の全関係と呼んではならない。

## 5. 保存量と未完了gate

一回のfollowupは0.614245秒、最大RSS111,608KiB。
176,910 scalar-DAG nodes、3,316,199bytesで60秒/2GB/32MiB内。
新しく行った数値演算は、保存したactual DAGの原点でのrational ADと48/9次の
小さいlinear algebraだけ。新しいbeta点・source点、GB、core/source再生成、series展開は無い。

専用source `code/ambient_B18_top_graph_telescoping_20261104_v63.py`。
保存 `reports/AMBIENT_B18_TOP_GRAPH_TELESCOPING_20261104_v63/` に、rank gate、全finite
inverse equations、全regular quotients、六relative graph pairs、十二top pairsを置く。
旧blocked reportは当時の正確なfrontierとして変更しない。

まだ未計算なのは、四module criterionの全lower係数と、この十二pairの実際のquotient moduleへの
所属、socle、depth19/20、full nilradicalである。係数の構成完了を、それらの解決と混同しない。
