# T37 の全固定 P1 からの正則 equivariant coisotropic planes

2026-11-08、PRIMARY `E_reduced_nonnormal`。純手還元。
既受理の整数実曲線 `K4A_SMOOTH_24_PARENT_ACTUAL_CURVE.md` の全十括弧を用いる。
追加の一般点、normal jet、period profile を選ばない。
以下の有限リストの actual section と磁気 defect は別の一括証明書で検査する。
配置だけから全六磁気源の帰属を推論しない。

## 固定 P1 と有限性

座標は `(c,a,b,u,v)`、torus weights は `(2,1,1,0,0)`。
P(0) の plane `span(c,Aa+Bb,Uu+Vv)` の coisotropy は正確に

```
U*[A*(U²+V²)-B*U*V]=0.
```

ここではこの全固定曲線の成分 `U=0` を扱う。
他成分は未処理として保持する。有限 chart の全固定 plane は `a=k*b,u=0`。
その全正則 equivariant lift は、projection が原点で同型であるため、正確に

```
a=k*b+A*e*c, u=B*e*b+C*e²*c.
```

基底の新 parameter を加えた形ではなく、固定した一つの actual curve 上の全 graph である。
e≠0 を torus で e=1 へ送れる。

`z=a-k*b-A*c,w=u-B*b-C*c` と置くと

```
{z,w}=Pau-B*Pab+C*Pca-k*Pbu+(A*B-k*C)*Pcb-A*Pcu.
```

plane への制限の bv,cv 係数は `6*k*(A-6B)` と `6*[A*(A-3B)-3kC]`。
k≠0 なら A=6B,C=6B²/k。K=k² と置いた残る全三係数は、非零定数を除き

```
2B²+2B+K*(1-15B),
6B²*(B-1)+K*(-45B²+33B-4),
6B⁴+K*(-45B³+63B²-26B+3).
```

K≠0 なので B≠0、15B-1=0 も第一式に反する。
K=2B(B+1)/(15B-1) を第二式へ戻すと `(5B-1)*(3B-1)=0`。
第三式は両根で零。従って k≠0 の全解は

```
A=6/5, B=1/5, C=k,  k²=6/25;
A=2,   B=1/3, C=3k, k²=2/9.
```

各二符号、正確に四 plane である。分母の補集合を捨てていない。
k=0 は既保存の全四 plane `(A,B,C)=(3,1,0),(2,2/3,0),(0,1,±sqrt6)`。

## chart at infinity

固定 P1 の残る点は b=u=0。この全正則 equivariant lift は

```
b=A*e*c, u=B*e*a+C*e²*c.
```

e=1 の `{b-Ac,u-Ba-Cc}` の av,cv 係数は `6*(6B+A),18C`。
従って A=-6B,C=0、続いて a² 係数135Bから B=0。
よって infinity には厳密に一つの lift `b=u=0` がある。
全 U=0 固定成分からの指定正則 equivariant class は、旧四、新四、infinity 一の全九 flag。

## actual normalcone と全磁気 defect の意味

各 plane の二 normal を z,w、base を三座標 y とする。
actual tensor で `{z,w}|S=0` なので、normal だけを e 倍し tensor を e 倍した
通常の normalcone pushforward は全負冪が消え、全 e の実 Poisson family になる。
その端点の mixed anchors は `P_{y,z}|S,P_{y,w}|S`、normal-normal は z,w に一次。
section は `det(E_y,P_{y,z}|S,P_{y,w}|S)`。

全 first-order plane 移動 `z=e*f(y),w=e*g(y)`、f,g の全六線形係数が作る磁気 defect は

```
f*∂z(Pzw)|S+g*∂w(Pzw)|S
 -sum_i(∂i f)*P_{yi,w}|S+sum_i(∂i g)*P_{yi,z}|S.
```

これは全 quadratic base 六係数への一つの線形 map である。
rank 六なら全磁気源へ直ちに到達するが、rank が落ちれば残商を勝手に零と置かない。
後者では親の真 period/edge 移動など別の実族が必要である。
section が triple factor を持つ flag だけなら、その全 finite class の陰性として止める。
全 T37、他の固定成分、singular base frame、非一様 weights は排除しない。
