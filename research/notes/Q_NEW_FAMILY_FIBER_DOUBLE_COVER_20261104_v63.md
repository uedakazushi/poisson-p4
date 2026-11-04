# 同じbeta上のactual有理fiberとJ直線の明示二重被覆

2026-11-04。明示Luroth生成元w=w_fiberを、新しいv=0 sectionの座標
(xi0,r0,s0,q0,B0)へ書き換える。旧t=2 sectionの同名座標とは同一視しない。
基礎体K_B=k(B_NEW15)、k=QQ,Cで、generic函数体の主張を行う。

```
h=s0²−25r0q0²,
A=200r0q0(s0−1), B=s0−25r0q0², D=r0q0²h,
J=(256/25)*(A+B*w)²/[D*(1600r0−w²)].
```

既存の五定数の境界値K=(8s0+u*q0)/(8s0−u*q0)、
Iv=u/[512(8s0+u*q0)]、u²=1600r0を、前packetのJ(w)へ代入するとこの式になる。
従ってactual曲線K_B(w)のJ直線への写像は次数2である。

関係式 J*D*(1600r0−w²)−(256/25)*(A+B*w)²=0 をwの二次式と見ると、判別式は

```
(80D)²*r0*J*(J−Jcrit),  Jcrit=256−256/(25r0q0²).
```

よって

```
K_B(w)=K_B(J)(sqrt(r0*J*(J−Jcrit))).
```

Jは超越的で、genericにJcrit≠0だから判別式は平方でなく、拡大次数は正確に2。
generic fiberの滑らかな射影模型は既知のP1で、その二つの分岐値は0,Jcritである。
分岐点もw=-A/B、w=-1600r0*B/AとしてK_B上に有理的である。
ここでgeneric baseからさらに除いた分母やJcrit=0の特殊層を、全域分類したとはしない。

実際のv=0境界sectionはw=200r0q0/s0を与え、上式によりJ=256となる。
この境界値だけからJ=256上の全点・全他成分への所属を推論しない。
NEWのbeta超曲面の定義方程式やspectral次数の主張でもない。

専用sourceは約分、u²=1600r0の剰余、二次判別式、分岐点と境界値をexact QQで検査する。
新点・rank・高次不変式・source26は使用しない。
