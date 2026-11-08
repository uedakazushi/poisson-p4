# 零symbol唯一軌道：有理Casimirと四つの二重零点

2026-11-08、PRIMARY E、plane_collision_review の純手計算。
入力は `PLANE_ZERO_SYMBOL_WHOLE_SINGLE_ORBIT.md` の既受理 canonical raw tensor のみ。
新 operator／kernel／profile／GB／jet の実行は零。既知親又は新成分性は未解決のままとする。

座標 (a,b,c,z,w) で全十括弧は

```
Pab=aw-4cz, Pac=-az+4bw, Pbc=bz+cw,
Paz=w², Paw=-z², Pbz=zw, Pbw=Pcz=Pzw=0, Pcw=-zw.
```

actual quintic は

```
F=3az²w²+b(z⁴-4zw³)+c(w⁴-4z³w).
```

## 1. D上の実際の有理第一積分

w=1、t=z での projective bracket は raw bracket ではなく、Euler補正を含む次式。

```
{a,b}=a+bt²-4ct, {a,c}=-2at+4b+ct², {b,c}=c,
{a,t}=1+t³, {b,t}=t, {c,t}=t².
```

t を反転する dense chart で q=b-c/t と置く。Leibniz ruleを全六式に直接使うと

```
{q,b}={q,c}={q,t}=0,   {q,a}=-F/t².
```

従って q は D上の有理Casimirである。全斉次表示は

```
q=(bz-cw)/(zw).
```

t=0 又は w=0 を新しい source complement として除外した命題ではない。
これは同じ有理関数／正則tensorの稠密開上の恒等式で、分母を払って全域へ戻る。
ambient rank四で q がCasimirとは主張しない。実際、上の最後の非零式がその相違を示す。

さらに Dを除く同じ chart で h=F/t⁴、r=-1/t とすると

```
{q,h}=-3h, {c,r}=1,
{q,c}={q,r}={h,c}={h,r}=0.
```

手確認に必要な中間式は {F,t}=0、{F,c}=-4tF である。
逆写像は

```
t=-1/r, b=q+c/t,
a=h*t²/3-q*(t²-4/t)/3+c*(t+1/t²).
```

これは有理Poisson座標表示であり、PGL同定又は既知閉成分所属ではない。
有理Darboux表示だけを actual quadratic family への変形へ読み替えない。

## 2. 全order二以上の零点はplane内の四点

F_a=3z²w²、F_b=z⁴-4zw³、F_c=w⁴-4z³w。
z,wが同時零でないとこの三微分は同時に零にならないため、被約Sing(D)は正確に plane S=(z=w=0)。
rank零点はSing(D)に入るので、plane外の新しい零点は無い。

a=1 chart で、plane上の点における projective bivector の一次jetは

```
{b,c}=(2b-4c²)z+(2c-4b²)w + terms of normal degree at least two
```

だけであり、他五括弧の一次jetは零。
例えば

```
{b,z}=2zw-bw²-4cz², {b,w}=w²+bz²-4czw,
{c,z}=-z²-cw²+4bzw, {c,w}=-2zw+cz²+4bw²,
{z,w}=z³+w³.
```

従って order二以上のprojective zeroは全斉次条件

```
ab=2c²,  ac=2b²
```

の全解。a=0では全b,cが零となり射影点でないので、a=1 chartが全解を覆う。
四点は

```
[1:0:0:0:0],
[1:ζ²/2:ζ/2:0:0]  (ζ³=1).
```

二quadricsのJacobianは各解でrank二なので四点はreduced。
三つのζ点のVandermonde determinant、及び原点と二つのζ点の2×2 determinantは非零で、
四点のうち三点共線はない。
本noteの帰属入口には「order二以上の四点」という弱い結論だけを使い、全局所同型型又はexact nilpotenceは主張しない。

これらは既T36/T38等の親幾何を照合する入口である。四点配置や有理Casimirが一致するだけでは
全150 tensorのactual PGL/Euler同定又はactual curveは得られない。
固定Cartan pencilの既失敗を再実行せず、同じ六法線のjet延長も行っていない。
