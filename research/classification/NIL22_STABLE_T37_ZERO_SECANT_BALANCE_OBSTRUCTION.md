# T37 zero-secant weighted balance：一つの全tensor pole障害と停止範囲

2026-11-08、PRIMARY E_multiple_hyperplane。純手検算、数値job無し。
T37のzero-secant geometryから選んだ一つのweighted classを全tensorで判定する。
source全period/smoothing積関係を保持する。全T37、別frame又は別period valuationの非所属は主張しない。

## 指定した実class

truecompleted T37のsourceを(a,b,c,u,v)、periodをp、smoothingをt0,t1,t2,t3とする。
zero-secant lineはa=u=v=0である。そのb,c非零部分の中心を対称frameに置き、
target順序を **(q,h,s,x,y)** として

```
q=epsilon²(b+c)/2, h=epsilon*a, s=(b-c)/2,
x=epsilon*u, y=epsilon*v,
tensor scale=epsilon.
```

inverseは

```
a=h/epsilon, b=q/epsilon²+s, c=q/epsilon²-s,
u=x/epsilon, v=y/epsilon.
```

GL determinantはこのtarget順序でepsilon^5/2。
同じ座標を(q,s,h,x,y)順に並べる場合は符号だけ反転する。
対象period classはpの非零定数項を保ち、t2,t3はregular、
`t0=alpha/epsilon+regular, t1=beta/epsilon+regular`、alpha,beta自由。
全regular tail、alpha又はbeta零も同じclassに含む。

一般secant点のD法線三次初項は、非零scalarを除き

```
a*(3p*b*u+(3/2)t2*c*a)*(2p*c*v-t3*b*a).
```

これは因子配置への入口だったが、Dだけでは実tensorの正則性は結論出来ない。

## 普遍的に消せないprojective係数

sourceの三つの必要entryは正確に

```
Pab=t0*u², Pac=t1*v²,
Pbc=p*b*c+t0*t1*t2*t3*u*v/(18p³)
             -t0*t1*t2²*t3²*a²/(48p^5).
```

後二correctionを省略しない。指定mapとtensor尺度を全て作用させたraw Πのentryは

```
Π_hs=(t0*x²-t1*y²)/2,
Π_qh=-epsilon²*(t0*x²+t1*y²)/2,
Π_qs=-p*q²/(2epsilon)+p*epsilon³*s²/2
      -epsilon*t0*t1*t2*t3*x*y/(36p³)
      +epsilon*t0*t1*t2²*t3²*h²/(96p^5).
```

従ってprojective sectionを表すE wedge Πの(q,h,s)成分

```
q*Π_hs-h*Π_qs+s*Π_qh
```

における **h*q²の係数は厳密にp/(2epsilon)** である。
Pab/Pacと元correctionはnormal変数だけから成るので、このmonomialを消す項は無い。
全smoothing積、全alpha/beta、全regular tailを保持しても同じ恒等式である。
p(0)!=0なので、このprojective poleは非零。Euler-vertical補正ではE wedge Πを変えられず、消去不能である。

この一係数は、指定tensor尺度での正則rank四curveが存在しないことの十分な全tensor障害である。
t0/t1の高次tailを足すことはこの係数へ影響しない。

## projective再尺度でもstable端点にならない

上のbalanced classでは、元十entryとinverseの次数から、全pushのpole次数は高々一。
負冪entryはs-row/columnだけに現れる。
実際、sourcePbcはpole高々四、mixed b/c-normalは高々三、normal-normalは高々二で、
mapとtensorを掛けるとq/h/x/y同士のentryは正のepsilon次数、sとのentryだけがpole一になり得る。

したがってepsilonをさらに掛けたprojective limit Kは

```
K=partial_s wedge W
```

というdecomposable raw bivectorである。
上のhq²係数によりE wedge Kは非零なのでprojective limitはrank二。
単にraw familyが非正則というだけでprojective極限全体を不存在と呼ばず、
その唯一の再尺度極限がstableのrank四端点と異なることまで区別した。

## 正確なscopeと停止

本classの全alpha/beta/tailは、両stable全rank四源の端点を作らない。
これは指定weighted zero-secant frameとp-unit valuationの結論だけである。
pを零へ動かす別valuation、異なるframe、別center又は全T37の非所属へ延長しない。
同classのperiod次数・tail・profileを増やす試みは停止する。
新matrix/operator、point、rank、minor、normal jet又はCASは実行していない。
次の義務はCartan/Borelの全source構造を保つ別のactual parent entranceである。
