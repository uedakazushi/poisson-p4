# p=0 冪零 cofactor の全五源は MPS T06 に入る

2026-11-08、PRIMARY `E_reduced_nonnormal`。全 coefficient producer PASS、独立 review 未完。
新しい成分は零。最小の陽性証明は下の literal 四 parameter profile と actual PGL だけであり、
探索に使った68未知数 matching、既 tangent 上界、旧 normal/jet を必要としない。

## 元の全源と endpoint

元の必要十分五源は `TWO_HYPERPLANES_CUBIC_SUPPORT.md` §5.2、p=0 の行。
自由変数 `n2,n3,n5,n8,n10` を全て残し

```
A=n10*n3+n5+n10²/2, B=n2+n3*n8,
U=u+A*h-B*x, V=v, W=w-n10*h-n8*x
```

とする。det=1、逆は `u=U-A*h+B*x,v=V,w=W+n10*h+n8*x`。
係数の局所化を一つも使わず、変更後は l=n3 として

```
{x,u}=xv, {x,v}=x(h+x), {x,w}=hx,
{u,v}=hu-v²/2+l*(-hv+(h+x)w),
{u,w}=hu-w²/2+l*h*(w-v),
{v,w}=h*(v-w)-l*hx.                                         (1)
```

全 l で cubic は `F0=(h(v-w)²+xw²-2hxu)/2`、actual D は `hxF0`。
F0 は u について一次、その係数 hx と残部が互いに素なので既約である。
全 source と normalizer/逆の150係数一致は producer で再生した。

## 真 T06 source からの literal curve

受理済み raw T06 の source 座標 `(x0,x1,x2,x3,x4)`、parameter `(p0,p1,t0,t1)` を使う。
全十括弧は `reports/MPS_S61_RAW_CATALOGUE_20261103_v62/T06.json` の actual family のままである。
parameter を次の四つへ同時に代入する。

```
p0=e+2e²l, p1=-e-4e²l, t0=e, t1=-e.                         (2)
```

実 PGL は

```
x0=-2u-2v/e-(h+x)/e²,
x1=x, x2=h, x3=w+h/e, x4=v-w+x/e, det=-2.                   (3)
```

逆は

```
h=x2, x=x1, w=x3-x2/e,
v=x3+x4-(x1+x2)/e,
u=-x0/2-(x3+x4)/e+(x1+x2)/(2e²).
```

e 非零で (2),(3) は真の T06 source 点への可逆 PGL である。
余分な tensor scale や coefficient product の独立化は行わない。
親の四 parameter は全値で projective Jacobi を満たす真の source で、
その全像閉包/PGL 飽和が既知の T06 component に入ることは既受理 catalogue を用いる。

三次の intrinsic な理由は、`C=UVW+UX²+VH²` が

```
X->X+tV, H->H+sU,
W->W-2tX-t²V-2sH-s²U
```

で保たれることである。(3) は t=s=1/e のこの actual automorphism を使い、
marked H,X を target h,x へ動かす。`2p0+p1=e` だから同じ cubic が全 e,l で残る。

## 正則 Euler 商と全境界

二次斉次 bivector P の標準代表を

```
C(P)=P-E∧div(P)/5
```

で取る。`div(E∧L)=5L-(tr L)E` なのでこれは Euler ambiguity だけを除き、
同じ projective bivector を表す。PGL は projective tensor として作用し、
source に (2),(3) を代入した C(P) の全150係数で e の負冪が消える。
e=0 は厳密に C((1)) と一致する。全 e,l について元の350 Jacobi も零。

実五次式は full signed Euler Pfaffian から直接

```
D(e,l)=(1+e*l)*(h+e*w)*(x+e*(v-w))*F0.                      (4)
```

従って e=0 は rank 四の (1) であり、l=0 と全有限 l が同じ正則 family に含まれる。
l や元の source 係数を一切反転しない。e 非零の一般点が既知 T06 に入り、T06 は閉なので、
全 l の endpoint は T06 に入る。元の全逆 normalizer により**全五源とその全 PGL 飽和閉包が T06 に含まれる**。
既知点の境界帰属を逆向きに使う推論や local uniqueness の仮定はない。

producer は全4 parameter の元350、全source normalizer/逆、(3) の行列逆と cubic 恒等式、
全 l endpoint150、全曲線350、実 Pfaffian を QQ 上で照合した（4.036秒）。
独立 reviewer は (2),(3) の最小 literal family だけを別算術で確認すればよく、
探索の68未知数線形系の rank や complete-profile 証明はこの陽性被覆の必要仮定ではない。

他の nilpotent cofactor、平面1+1+3の残三配置、Ga9 linear-w、stable二源などの全 E 残余は別に保持する。
