# T/D12の補集合：八tailの同一source族

2026-11-05。SOURCEと支持・高さ・生存期間の前見積だけであり、新しい列の値やκ内積を
計算しない。点・標数・GCD・factorization・rank・minor・旧572sourceの再実行も行わない。

現時点で独立受理済みの新しい十分開は、FULLUの`D(T*D12)`である。
TはKAPPA_EXTRA2_ENTRYのroot239、7208項、a指数最小4。従って既知単元a四乗を
取り除くことはできるが、C等の未知の割切りやTとD12の共通因子を仮定しない。
既知恒等式`kappa3*D8=cprod*Saux*T`は因子の積の恒等式であり、T=0を既に処理した
Theta/Bbad集合と同一視する根拠にはならない。

## 既存の全点被覆との比較

受理済みcoverageは`750000p=g1*Theta`と
`1=Theta/B²+C*g2*g1/B²`でD(p)、g1=0、Theta=0を扱う。
actual rawH rank<7の幾何学的支持は、actual Bbad13とTheta上のP80/P100二曲線。
それぞれでfull cubic57が独立に証明済みである。特にheavy rankの欠損からfull cubic欠損を
推論してはいけない。四次239の結果やP100の有限例外テストもここで全曲線へ拡張しない。

Q/Lの二unit completionは全FULLUのsource変換であり、単独では新しい不良集合の消去を
意味しない。古いG0=Q*L四乗の除去はactual Bbad13スキーム同定に含まれる依存であって、
任意のQ=0/L=0点でfull57という別定理ではない。
したがってT=0/D12=0に残り得る実欠損はH7開の中に限られるが、既存三集合だけでは
これら補集合を消したことにならない。

## 同じM12を保つ八つのdirect15恒等式

current15の固定四行は`[2,4,5,8]`。H2列`[27,31]`とanchor22の4×3行列Aについて、
`kappa_i=(-1)^i det(A_without_row_i)`という**既計算の四多項式**を使う。
anchorはliteral M[:,22]であり、余分なaを掛けない。frameは受理済みtail2計算そのもの。

九pure列`[22,11,13,25,15,19,14,18,17]`上の四行のrank一以下と、
H3三列`[28,30,29]`のliteral零は既に受理済みである。新しい八tailは
`j=[1,2,3,4,6,7,8,9]`。各列のa-Laurent/QQ/B分母をsourceから型付けし、
四行で共通なB冪b_jだけを払って

```
T_j=B^b_j * kappa * M[:,j]
det([A,B^b_j M[:,j]])=-T_j
```

と置く。T_2は既存Tである。
M12は全く変えず、DIRECT_M12の十五角の追加tail2だけをjへ替える。
source行順とH2/H3/L/gamma/aの倍率を維持し、j列だけB^b_j倍する。
既知九pure列にjは入らないため、どのjでも十五列は重複しない。

三つのpair行をanchor係数で消す同じ多項式Laplace証明により

```
A_anchor * D15_j = -T_j * D12
```

となる。Hcorrのdetは`A_anchor²*det(A4_j)`であり、両辺のA二乗をsource domainで消す。
Aやκ3を追加chart単元にしない。従って各jで`D(D12*T_j)`はfull57の十分開である。
八つ合わせると

```
D(D12) intersect union_j D(T_j)
```

を得る。このsource族はT=0を扱うための具体的な新しい計算目標である。
T_jたちがunit idealを生成する、実bad Fitting idealと等しい、あるいは九pure列以外の
rankが分かるとは主張しない。全T_jが零であっても、元のphiで失ったpure方向により
actual full57が回復する可能性がある。

## 実行候補の選び方と資源

全八列についてactual Q2/L finite DAGのsource閉包、既知a/B分母、exact QQ分母、
指数box・total degree・整数l1上界だけを伝播する。新しい多項式積は作らない。
κ四個は保存済み多項式を再利用する。全四行と全八tailの効果を維持する。

将来の一列計算には、四つの積を別々に実体化しない単一整数accumulatorを提案する。
source列をlast-useで作ってから、共通QQ分母の下で四内積を直接加算する。
coefficientとpacked exponent keyの双方、一原子分のJSON数値metadata、15000 work、
source読込、出力・readbackを含めて見積る。scriptは最小候補と全八profileを保存するが、
これ自体はworker実行の承認ではない。実workerのcustody・guard・source inverseは別gateとする。
選択されたj以外への自動拡張、別minorや新しい点を使う代替はない。
