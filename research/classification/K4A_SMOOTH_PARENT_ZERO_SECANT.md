# K4a actual24 親の零直線と二つの distinguished points

2026-11-08、PRIMARY E、印刷された整数Pi(1)からの純手計算。
T37標準零曲線のframeに使う具体データであり、親のPGL同定自体はまだ主張しない。
新zero operator、point列又はprofileを実行しない。

## Actual rank-zero line

`K4A_SMOOTH_24_PARENT_ACTUAL_CURVE.md` のe=1のtensorは

```
ell: a=3c, u=b, v=9c                                       (1)
```

上でprojectiveに零となる。これはcoisotropic plane `(A,B,C)=(3,1,0)` のsection
`c²(v-9c)³` の三重直線である。
そのplaneでv=9cも入れると全十entryは

```
Pca=-60bc, Pcb=-2b²-801c², Pcu=-2b²+9c², Pcv=0,
Pab=54b²-2403c², Pau=54b²+27c², Pav=540bc,
Pbu=810bc, Pbv=18b²+7209c², Puv=18b²-81c².
```

全六entry `T_cij=c*Pij-x_i*Pcj+x_j*Pci` が恒等零。
例えばT_cabは `(54+6-60)cb²+(-2403+2403)c³=0`、
T_cbvは `c(18b²+7209c²)+9c(-2b²-801c²)=0`。
c非零なら全T=E∧Pが零、閉性からc零の点も含めてell全体がrank-zeroである。

## 全直線上のfirstjet rank-drop

affine c=1、line coordinate b、normalを

```
z=a-3, w=u-b, h=v-9
```

とする。actual affine Poisson matrixは `pi_ij=Pij-x_i*Pcj+x_j*Pci`。
六entryのfirstjetのうち三つを pi_ab,pi_bu,pi_bv に取り、
他三つは同じrow又は零になる。normal方向の全firstjet matrixは

```
       z                w          h
[ 6b²-153,              0,        54 ]
[ -108b,        -9-6b²,        36b ]
[ 432-36b²,             0,   6b²-153 ].                     (2)
```

determinantは `-27*(2b²+3)³`。
従ってc非零のell上でfirstjet rankが3より落ちるのは正確に

```
q_±=[c:a:b:u:v]=[1:3:±sqrt(-3/2):±sqrt(-3/2):9].             (3)
```

そこでrankは一である。
無限遠q_inf=[0:0:1:1:0]はaffine b=1へ移す。
normal `(a-3c,u-1,v-9c)` に対しpi_ca,pi_cu,pi_cvの三行は
`(6,0,0),(0,6,0),(36,0,-6)`、det=-216だから追加rank-dropは無い。
これで直線上の全例外点を保持する。

q_v=[0:0:0:0:1]もactual rank-zero点でellの外にある。
T37の標準quartic+secant frameを作る際、ell、q_±、q_vは非任意の実幾何データとなる。
ただしquartic成分そのものの全parameterization、T37への全150inverse、
又は全親の帰属は別のactual等式を要し、本ノートだけでは結論しない。
