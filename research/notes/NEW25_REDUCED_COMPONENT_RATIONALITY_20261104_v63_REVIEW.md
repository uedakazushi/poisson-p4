# NEW射影被約成分25の有理性：合成証明の独立監査

2026-11-04。対象は既存の幾何的既約成分 `R4_POSITIVE_WEIGHT6_25` の被約台。
full-tensor fiber、実S61移送、B_NEW15有理性の三つのFINAL独立gateと、
それらを適用する同名author packetのhashを専用reportで照合してから結論を受理する。
新しいfull tensor展開、source jet、不変式、点、prime、GBは計算しない。

## 三つの異なる空間と同じbase

`B=B_NEW15` はprimitive firstjetの被約像、`Y=Y16` はactual full S61の被約成分、
`X=X25` は元Poisson空間の射影被約成分である。
firstjetだけの有理性からXの有理性を推論してはいない。
受理済みB_NEW15応用gateは、actual invariant pullbackの指数1、有理section、
全幾何的stabilizerの自明性を使い、QQ上で
`B ⇢ Sp4×Spec E5`、`E5=QQ(xi,r,s,q,B)` を双有理に同定している。
従って `K=QQ(B)` は純有理十五変数体である。

raw pairの整合性にも独立な根拠がある。受理済みraw kernelのperspective成分が恒等零なので、
source6+GL4のraw像はdim21、値Jへの正規化後のpair像はdim15。
t=2の有理pair sectionのSp4飽和はB15へdenseに射影し、同じdim15のpair像にもdenseである。
`Sp4×S5→pair像→B15` の合成が既受理の双有理写像なので、中間のpair像→B15も双有理。
従って一般部のβ matchingはαも合わせ、βのinvariantだけからαを推測していない。

full-tensor gateは、raw pairのorbitを合わせる有限有理frameを作り、
そのframeを**全quadratic tensor**へ作用させる。primitive firstjetのkernelを
full tensorの対称性とは仮定しない。t依存のquadratic項を残す。
さらに既知の有限tau同値はfull tensorに対して成立し、
`F6=E5(u,t)`、`u²=1600r`、`tau(u)=-u`、`tau(t)=t/(t-1)` のもとで

```
z=u(t-2)/t,   t=2u/(u-z),   tau(z)=z,
F6^tau=E5(z)
```

を与える。uをQQ(E5)内の有理函数としていない。
最後の固定体表示は、F6=E5(z,u)でtauが二次拡大の非自明な符号作用となることによる。
raw pairのmatchingだけでfull tau固定を代用しない点が本質である。

得られる支配的写像は `Sp4×Spec E5(z) ⇢ Y` であり、firstjetへの合成は
**同じ** `Sp4×Spec E5 ⇢ B` の双有理写像とcompatibleである。
physical pushforwardのgroup変数gに対するβはrho(g⁻¹)なので、h=g⁻¹という
有理なgroup inversionで既受理のpullback規約rho(h)へ揃える。
よってpullbackは単なるQQ体の包含ではなく、同じK上で

```
K ⊂ QQ(Y) ⊂ K(z)
```

となる。異なる有限spectral base上の有理被覆からの結論ではない。
元sourceのfull像とSp4飽和がactual Yに稠密であることもfull-tensor gateに含まれる。
parameter依存の可逆なprojective操作が同じPoisson成分を保ち、全full像の次元16を
保持することがここに必要である。

## LürothとQQ降下

dimY=16、dimB=15なのでQQ(Y)/Kは超越次数1で、Kへの定数写像ではない。
任意の体K上のLüroth定理は、`K⊊L⊂K(z)` に対し `L=K(z0)` を与える。
従って `QQ(Y)=K(z0)`、YはQQ上で有理16次元である。
zがその生成元であることや、元の一変数写像が次数1であることは必要ない。
z0の展開を今回構成したとは主張しない。

全dominant mapとfield inclusionはQQ上で構成されているため、C上の有理性を
QQへ逆降下していない。幾何的既約性は元の幾何的source像とaccepted component gateにより
保持される。QQ上の双有理写像をCへbase extensionしても双有理であり、C有理性も従う。
generic曲線はこの関数体の意味で幾何的integralだが、各特殊fiberや全Jacobi schemeの
被約性を意味しない。

## 実S61から元射影成分へ

移送gateは固定評価点での非退化値Cを六つの係数として記録し、有理Darboux行列を適用、
実旧S61の四つのperspective gaugeを除く。逆は記録した四係数とCを用いる有限変換である。
accepted55消去と134式のideal恒等式により、形式近傍だけでないactual scheme開の

```
U6×A4×S61
```

を得る。平方根や余分な有限被覆を挿入していない。
幾何的連結PGL5とscalar群は各幾何的既約成分を保つ。
よって両方向の操作は対象成分を保ち、積表示をその被約成分とactual Yへ制限できる。
QQで定義された幾何的source閉包なので、この制限もQQ上で成立する。

projectivizationは正規化**前**に値の一係数c=1でscaleを固定する。
固定評価点を動かす四自由度を加えることも、S61上でscalar作用を再び商にすることもない。
従って

```
X ⇢ U5×A4×Y,
QQ(X)=QQ(Y)(t1,...,t9)=QQ(t1,...,t25)
```

となる。Y16への追加自由度は正確に9である。

## 受理範囲

結論は既存NEW射影**被約**成分のQQ/C有理25次元性であり、成分数を増やさない。
generic full Poisson完成環の長さ5とCM type2を変更しない。
全Jacobi schemeが被約または滑らかだとは言わない。
NEW firstjet像の実定義式、正確な次数、全boundary軌道、一般betaの五点尽くし、
全Poisson成分の上界・完全分類は未解決のまま残す。

専用checkerはFINAL gateとauthor packetの保存hash、範囲flags、次元算術を照合する。
author producerの再実行や新しい算術探索はしない。
