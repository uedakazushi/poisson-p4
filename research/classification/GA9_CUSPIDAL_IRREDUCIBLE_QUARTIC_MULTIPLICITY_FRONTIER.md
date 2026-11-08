# Ga9 cuspidal24の既約quartic境界ではmarked点の重複度は二

2026-11-08、PRIMARY `E_reduced_nonnormal`。DVR上の三係数恒等式だけによる補題。
新CAS、点又は局所jetを用いない。Ga9 linear-w全源の親候補を一つ除外する。
本補題の独立reviewは未完であり、既知成分及び元の全sourceの受理済み状態とは区別する。

## 1. 既知の全源からproperなmarked flagへ

C=`R4_GA9_CUSPIDAL_PLANE_QUARTIC_24`とする。
その全源と成分性は `GA9_CUSPIDAL_QUARTIC_COMPONENT.md` 及び独立reviewの受理済み入力である。
一般の実退化式はD=hF_aで

```
F_a=(hW+6y²)²−144y³(h+ay).
```

H={h=0}、p=[0:0:1]とするとF_a|Hは非零y⁴の倍で、pのquadratic coneはh²。
flagに適合する一般quarticを

```
F=a0 h²W²+hW B2(h,y)+C4(h,y)
```

と書く。generic sourceではa0≠0で、binary quartic

```
G4=B2²−4a0 C4=576y³(h+ay)
```

は三重根linear formを持つ。
H、F、p及びこのlinear formのprojective flagはproperである。
exact essential3開で商平面を正則に固定すれば、Cの境界への任意のspecializationにこれらのflagも延びる。
flag座標のW→dW+L(h,y)でG4はd²倍されるので、三重根の条件は座標に依存しない。
sourceの稠密開からのcurveを取り、必要なら有限DVR拡大して三重根linear formを選べる。

## 2. 既約quarticのa0=0境界は実際にはliftしない

targetでもFが既約で、F|Hが非零の四重接触を持つとする。
a0=0を仮定して矛盾を導く。
B2=0ならF=C4はbinary quarticなので複素数体上で可約。従ってB2≠0。
target G4=B2²が三重根を持つことから、B2は一つのlinear form lの平方の非零倍である。
もしC4|l=0ならhWB2+C4全体がlで割れ、Fの既約性に反する。
したがってC4|l≠0である。

generic a0≠0のsourceからtargetへspecializeするDVR Rを選び、三重根lをproper flagで延長する。
lと独立な補linear formをxと書く。R上の係数を

```
B2=A l²+B lx+C x²,
C4=d0 x⁴+d1 lx³+d2 l²x²+d3 l³x+d4 l⁴
```

と置く。targetでA及びd0は単元、B,Cはmaximal idealに属する。
三重根条件G4∈(l³)のx⁴,lx³,l²x²係数は厳密に

```
C²=4a0 d0,
2BC=4a0 d1,
B²+2AC=4a0 d2.
```

generic a0≠0とd0単元からCはRの非零元である。
最初の二式よりB=C*d1/(2d0)。第三式へ代入し、integral domain R内で非零元Cをcancelすると

```
2A=C*(d2/d0−d1²/(4d0²)).
```

targetで右辺は零、左辺は単元となり矛盾。
Cを単元としてtargetを捨てたのではなく、R内の恒等式の非零因子をcancelして得た等式をspecializeした。
唯一の係数反転d0はtargetで非零と既約性から正当化されている。

従ってCの**単純H+既約quartic**のこのmarked境界ではa0≠0であり、pの重複度は正確に二である。
単にclosed三重根方程式を代入するとa0=0の余分なvertical解が存在するが、
上のDVR議論はそれらがgeneric sourceの閉包から来ないことを示している。

## 3. Ga9 c0 linear-w全源への適用

`HYPERPLANE_QUARTIC_NIL211_SHORT_TRANSLATION.md` §8及び独立
`NIL211_SHORT_TRANSLATION_REVIEW.md` の完全linear-w sourceは

```
F=k h⁴−L h³W−D4(h,y),
L≠0, D4|h=0=(r/4)y⁴, r≠0
```

の実quarticを持つ。
これは全parameterでprimitiveなW一次式で絶対既約、本質変数数は正確に3である。
唯一のHはh=0、四重接触点はp=[0:0:1]。
そのquadratic coneは零で、a0=0、B2=−Lh²≠0、C4|h=0=−(r/4)y⁴≠0。
よって§2の矛盾が全parameterで適用できる。

従って**Ga9 linear-wの全十一parameter源と全PGL飽和は、Ga9 cuspidal plane-quartic24成分の外**にある。
六磁気係数はこのFを変えず、磁気零も含む。
r=0又はL=0は元の既約quartic sourceの補集合であり、そこに本除外を無断で適用しない。

既知の全global T39/T40への非所属とは別のwhole-divisor幾何結果である。
T36/T37/T38、未知ambient成分の存在・次元、又はGa9 linear-wの全E帰属は未解決のまま保持する。
