# NEW25 source kernel の五 first integral と定数体の独立証明

2026-11-03 の継続作業。対象は六変数 **source の微分体**である。実際の NEW
beta image の Sp4 不変関数体と同一だとは主張しない。有限指数の問題は残る。

受理済み30×26 kernel の source 部分を K_raw とし、非零有理倍
`delta=(v/3072)K_raw` を取る。定数体はこの変更で変わらない。

```
r=384vX+256+v(3072W−1), s=384vX+256−v(3072W−1)
X=(r+s−512)/(768v), W=(1+(r−s)/(2v))/3072
```

元の全 polynomial kernel から次の6恒等式を独立に検査する。

```
delta u=0
delta r=r(r²−256²), delta s=s(s²−256²)
delta v=v(r²+s²)+16u(r−s)²
delta c=c(2*256²−r²−s²−32u(r−s)²/v)
delta B=B(2*256²−r²−s²−48u(r−s)²/v)
```

z=1−256²/r²、z_s=1−256²/s² とおく。両者の logarithmic derivative は
`2*256²` である。従って以下の五つの導関数は全て厳密に零である。

```
u,
k=z/z_s,
Iv=[v−8u(r−s)²/(rs)]/(rs*z),
Ic=c*v²/(rs*z²),
IB=B*v³/((rs)²*z³).
```

代数的独立性と定数体全体の主張も単なる次元数からは出さない。
`256²(1/r²−k/s²)=1−k` は k 上の有理点 `(1/r,1/s)=(1/256,1/256)` を持つ
conic である。その点からの line parameter
`t=(1/r−1/256)/(1/s−1/256)` の逆は

```
r=−256(t²−k)/(t²−2kt+k),
s= 256(t²−k)/(t²−2t+k).
```

両向きの有理恒等式を検査し、さらに

```
v=Iv*rs*z+8u(r−s)²/(rs),
c=Ic*rs*z²/v², B=IB*(rs)²*z³/v³
```

から全 source の逆を得る。従って五個の first integral と t は代数的独立で
`QQ(u,v,c,X,W,B)=QQ(u,k,Iv,Ic,IB)(t)` である。また
`delta t=256(r−s)t` は非零。この体では delta は定数係数体を殺す非零な
`f(t)d/dt` なので、標数零でその定数体は正確に `QQ(u,k,Iv,Ic,IB)` となる。

著者の別の conic parameter `(r/s−1)/(r/256−1)` とその両向きの逆も検査する。
その t=2 section は `v=−98304X/(147456X²−(3072W−1)²)`。これは五つの変数
を自由に残す記号的な section で、数値点の評価ではない。受理済み六変数
回路の全 node をこの置換で写し、保存された1032-node回路の全出力と照合する。

さらに既知の involution `(-u,-v,-c,-X,W,B)` は、raw value と全一次係数の
双方において `diag(-1,1,-1,1)` による実際の GL4 同値であることを検査する。
first integral 上の作用、section の保存も確認する。この section の関数体
`QQ(u,c,X,W,B)` の既知 involution による固定体は正確に
`QQ(X²,u/X,c/X,W,B)`。大きい体はこの固定体に X を付加した二次拡大で、
非自明な involution が X を反転するためである。これは**既知の involution**
だけの固定体であり、実際の全 orbit 同値がこれで尽くされたとは言わない。

この結論は有理微分体の一般点のもの。表示の分母が消える境界への一律の
延長、全 source fiber の同値性、完備な Ga 作用、実際の beta invariant
quotient への field index 1 は含まない。新しい beta 点や高次数 probe は無い。

```sh
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_new_family_source_constant_field_review_20261103_v63.py --output work/source-constant-field-independent
```

60秒・2.2GB上限。これは元の kernel に直接束縛した独立恒等式検査である。
