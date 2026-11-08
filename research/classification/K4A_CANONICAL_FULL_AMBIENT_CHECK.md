# K4a canonical 全一軌道：一つの全ambient成分性判定

2026-11-08、PRIMARY G/E。全sourceの canonical orbit は
`PLANE_CONCURRENT_M3_K4A_ORE_NONDOMINANCE.md` §§1–5 に純証明されたもの。
その whole-module theorem の独立受理と、以下の一つの canonical point の成分性を分離する。
原始同じ点の全350微分と全51作用列だけを一回取る。D≤101 なら接上界と障害を保存し、
追加 normal/jet/別point に進まない。

## Canonical actual tensor

座標は `(c,a,b,u,v)`。`P=Q+∂c∧W`、

```
Qab=av, Qau=u², Qav=2uv, Qbu=uv, Qbv=u²+2v², Quv=0,
Wa=2ab/9, Wb=3a²/2-b²/9,
Wu=-av/3+2bu/9, Wv=-4bv/9+4au/3.
```

計算では P を18倍した同じ projective point の原始整数代表を使う。
actual D は計算から正確に取り、以下の primitive polynomial の非零定数倍と照合する。

```
F=u⁴c+a²v³-2abuv²+(b²+3a²/2)u²v-4abu³/3.                  (1)
```

u⁴係数と残部は互いに素、cについてprimitive linear だから C 上絶対既約。
`∂cF=u⁴` より被約 singular support は u=0 上にあり、v非零なら `∂aF=2av³` からa=0。
逆に全偏微分は `{u=a=0}` と `{u=v=0}` の二平面で零なので、この二面で尽くされる。
両者はline `{u=a=v=0}` で交わる。

## 最弱の成分性矢印

全150二次係数に対する元350 Jacobiの微分 D と、
Euler25/GL5作用25/scalar1 の実51列 M を定義する。
Eulerのrankは24、D*M=0を全QQ係数で確認する。
各行列について全QQ rank、独立行列の非零determinant、全疎行列を保持する。

rankD=102、rankM=48なら

```
projective tangent dim=150-102-24-1=23,
actual orbit dim=48-24-1=23.
```

従って同じ点はsmooth unique23、actual orbit closure が全ambient成分である。
この場合だけ全source one-orbit theorem と結んでwholeK4aをその成分へ覆える。
既知126との相違と別作者 ambient 再構成が完了するまで新成分は計数しない。

`(a,b)->s(a,b),(u,v)->t(u,v),c->s²c/t` と tensor scale `t/s` は元Pを保つ。
射影で一dimtorusなので、全orbit次元の上界23が事前にある。
この torus だけで正確な orbit 次元又は成分性を推論しない。

## Ga9既約五次との intrinsic な差の手補題

二singular planesの共通lineの中に direction e を取り polar `∂eF` が非零四乗となるかは、
projective hypersurface の不変条件である。K4aでは共通lineの c direction が u⁴ を与える。
既知 Ga9 distinct quintic の二面の共通lineは C,E directions、polar は

```
α*(160A²BD-B⁴)-240β*A²B².
```

これが非零linear-form四乗なら、Bが割るのでそのlinear formはBに比例する。
A²BD係数からα=0、次にA²B²係数からβ=0で非零に反する。
従って両 hypersurface は projective 同型でない。
nil22-long/K22 の三 singular planes とも面数で異なる。
他の既成分の全相違監査は separate accepted-source factor/essential-variable 証明を用いる別義務である。
