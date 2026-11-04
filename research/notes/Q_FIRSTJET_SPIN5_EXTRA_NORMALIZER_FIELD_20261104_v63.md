# 追加 Spin5 成分の六座標作用と固定体の有限回路

2026-11-04。受理済みの連結 Spin4 joint biform 不変式体 K=C(I0,…,I5) に、
元の Spin5 normal-line stabilizer の**別成分**を作用させる。
旧 U0-line normalizer の μ4 とは異なる。純有理性は未解決であるが、
追加 involution の六つの像と固定体の13生成元を有限な有理回路として与える。
β 点、高次数不変式点評価、GB は使わない。

## 1. 追加 involution と交代する δ

既受理の元12係数の作用は

```
τ:(U(x,y),V(x,y)) ↦ (−V(y,x),U(y,x)).
```

二乗は (−U,−V) で、連結 Spin4 の中心の作用である。従って K 上では τ²=1。
cross4×4行列の determinant δ は τδ=−δ、δ²=16p0 を満たす。
保存された元の20項 δ と12×12作用行列から符号を全係数で確認する。

旧受理座標の z=c² と書き、v1=1 の torus gauge を取ると

```
r=zI1, p=(I2+zI3)/2, q=(I2−zI3)/(2r³),
s=(I4+zI5)/2, t=−(I4−zI5)/(2r⁵), z²=I0,
V=x0(p y0²+y0y1+t y1²)+x1(s y0²+r y0y1+q y1²).
```

ここで r,p,q,s,t は旧 torus 不変座標の gauge 表示である。
R=Res(V) は

```
R=(pq−st)²−(pr−s)(q−rt)
```

で、z²=I0 により六つの I の有理式になる。全式を action 証明書に保存する。
また

```
δ=[I0³ I1⁶+I4²−I0 I5²]/(4I0² I1⁵)
```

を厳密に得る。特に δ は K 内の非零 anti-invariant。

## 2. 交換後を正規化する exact frame

交換後の U' は −V の二つの quadratic row、V' は c(x0y0²+x1y1²)。
J0=2(pr−s)、J1=4(pq−st)、κ²=R と置く。一般係数体では J0,R は非零。
次の A は determinant1である。

```
A00=1/2−J1 κ/(8R),  A01=(−J1−4κ)/(2J0),
A10=J0 κ/(4R),       A11=1.
```

U'(Ax,y) の二つの quadratic は mixed x0x1 係数を持たない。
その pure-square 係数行列を M とすると detM=−κ。
この四つの恒等式は、p,q,r,s,t を全て独立とした QQ(p,q,r,s,t)[κ]/(κ²−R) で
直接確認する。第四根を有理関数と扱っていない。

z'=c'²=−κ とし B=c'M⁻¹ とすれば detB=1、U'(Ax,By)=c'U0。
N=M⁻¹ と書き、W の六係数を

```
W=(A00x0+A01x1)(N00y0+N01y1)²
 +(A10x0+x1)(N10y0+N11y1)²
 =x0(w0y0²+w1y0y1+w2y1²)+x1(w3y0²+w4y0y1+w5y1²)
```

で定めると、交換後の V の係数は c z' w_i である。従って六つの像は

```
τI0=R,
τI1=−zκ w1w4,
τI2=I0 R²(w0w1³+w5w4³),
τI3=−I0 Rκ(w0w1³−w5w4³),
τI4=zI0 R³(w3w1⁵−w2w4⁵),
τI5=−zI0 R²κ(w3w1⁵+w2w4⁵).
```

## 3. これは六変数上の有理回路である

上の計算を rank4 の有限 étale K-algebra

```
E=K[z,κ]/(z²−I0,κ²−R)
```

で行う。I0,R は非零に局所化する。二つの quadratic 拡大が独立な体になるとは
仮定しない。split や従属の場合にも rank4 étale algebra の trace は使える。

κ の符号は交換後 U' の二つの critical root の順序を変え、その正規化の二選択は
受理済みの**全** U0-line normalizer で結ばれる。六つの I はその全群で不変なので
像は κ の符号に依らない。z の符号は旧 cross-section の正規化選択を変えるが、
τ は実際の元表現の作用であり、その後の全 normalizer 不変量は同じである。
従って六つの像は E の K 部分に属する。

basis(1,z,κ,zκ) の係数1を取ることは trace_E/K/4 に等しい。
producer は各係数を**六変数上の四則演算 DAG**として保持し、六つの像について
この係数を出力する。これは単なる algebraic root 回路ではなく、base K 上の
有限有理回路である。三つの非定数 basis 成分を巨大な有理式へ展開して零と
比較したとは主張しない。その消滅には上の全 root-choice 不変性を使う。

τ²=1 はこの回路の巨大な自己合成から推測するのでなく、元12係数での二乗が
連結群の中心に等しいことと、回路が実際の作用を計算することから従う。

## 4. 固定体を生成する13有理式

F=K^τ とし、六つの Xi=Ii について

```
d=δ²,  Si=(Xi+τXi)/2,  Ai=(Xi−τXi)/(2δ)  (i=0,…,5)
```

を取る。全13式は固定され、Xi=Si+δAi。従ってこれらの生成体 F0 に対し
K=F0(δ)、δ²=d∈F0。τδ=−δ≠δ だから [K:F0]=2、ゆえに F0=F。
全13式を同じ DAG の出力に保存する。

13生成元は代数的独立とは主張しない。F の超越次数は6だが、六つの純超越基底へ
有理的に削減できるかは未解決である。特に、rational field の involution の
固定体が自動的に rational になるという推論は使わない。

既受理の associated-bundle 同定により、W16 の p0≠0 上の Sp4 不変式体は
この F に一致する。関数体の同一性なのでこの開への制限で十分である。
多項式不変式環、全 quotient の regular presentation、NEW25 の式は求めていない。

## 5. generic freeness と stable rationality の上界

generic freeness は別に正確に確認できる。p0≠0 では maximal minors の共変ベクトル
n が非零・非等方である。元 g が Ψ を固定すれば、lineだけでなく n 自身を固定する。
従って stabilizer は connected Spin4 に入る。一般 U=cU0 の exact stabilizer は T、
一般 V の weight1 係数 v1≠0 は t=1を強制する。これらの開条件は独立係数上の
非零条件であり、δ=c²(v1v4/4−v2v3) も恒等零でない。よって一般 Sp4 stabilizer は自明。

したがって W16 の rational quotient の一般 fiber は Sp4-torsor。
任意の標数0の体上で4次元非退化交代形式は symplectic Gram–Schmidt によって
標準形になるため H¹(k,Sp4)=1。その torsor は自明で、
W16 は quotient×Sp4 に双有理である。Sp4 の Bruhat cell は4+2+4=10次元の
有理開集合。ゆえに **F に10個の独立変数を付ければ純超越体になる**。
これは stable rationality の level 上界10であり、純有理性や最小levelを意味しない。

## 6. 制限した失敗と再現

全五像を QQ(p,q,r,s,t) 上で同時に展開する一度の試みは60秒で停止した。
frameと六つの W 係数までは0.316秒で完了していた。停止状態とその source を
report folder に保持し、同じ展開を再試行していない。
新 producer は有限 DAG と上記の小さい exact frame 証明だけを作る。

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_firstjet_spin5_extra_normalizer_field_20261104_v63.py
```

残る課題は、13式間の関係を使って F の純有理性を証明または別の障害を示すこと、
および NEW25 がこの6次元商のどの divisor に対応するかを同定することである。
