# Stable 二全源の共通 smooth quadric zero surface と全 parameter の実正規化

2026-11-08、PRIMARY E、既全 source からの純手恒等式。
親 component の認識、全 zero scheme の分類、又は新componentを主張しない。

座標を `(h,q1,q2,x,y)` とする。既全六磁気源の双方で、h=0への raw tensor restriction は

```
Pq1,q2=R2=m3*x²+m4*x*y+m5*y²,
Pq1,x=-x²/3, Pq1,y=-xy/3,
Pq2,x=-xy/3, Pq2,y=-y²/3, Px,y=0.
```

line_conic の追加 h(q1+q2) 項及び全低位磁気係数はこの restriction では零。
E_H=(q1,q2,x,y)、V=x∂q1+y∂q2 と置けば、全六 entry の恒等式は

```
P|H = (E_H∧V)/3 + (Q/3)*∂q1∧∂q2,
Q=q2*x-q1*y+3R2.
```

従って E_H∧P の全四三次 entry は、非零定数を除き `(xQ,yQ,0,0)`。
h=0内の reduced projective rankzero support は正確に smooth quadric surface S={Q=0}。
Qの q/base cross block は非退化なので、任意の m3,m4,m5 で rank四quadricである。
zero support が vertex lineだけだと解釈してはならない。

全 source parameter 上で actual linear shear

```
q1'=q1-3m5*y,
q2'=q2+3m3*x+3m4*y,
(h',x',y')=(h,x,y)
```

は determinant一、inverseは opposite shearであり、Q=q2'*x'-q1'*y' に戻す。
一つも parameter を反転せず、全 magnetic zero complement を含める。
これは zero surface の全実正規化であって、tensor の全 residual magnetic class を零へ送ったという結論ではない。
同じ shear は低位の磁気係数を変え得る。

scheme restrictionにも短い記述がある。Hのcoordinate ringで全zero idealは

```
(xQ,yQ)=Q*(x,y).
```

Q∈(x,y) なので S 上の vertex line ell={x=y=0} に embedded correction があり、

```
0 -> O_ell(-2) -> O_(Z restricted to H) -> O_S -> 0.
```

これは H 内の restriction の記述であり、h方向を含む full projective zero scheme の同定ではない。
双方の全源が同じ S/ell 配置を持っても、normal h derivative と residual magnetic quotient は異なり得る。

Euler ambiguity を除けば H 上の constant skew係数は
`A=(1/3)*∂q1∧∂q2`。A は H-tangentかつ image plane が Q-isotropicである。
従って既 quadric-cone source の H-normal coefficient 非零による既約D正規化を適用出来ない。
ここは H-normal coefficient 零から H因子を持つ、まさに可約Dの補集合に当たる。
同じ zero surface、又はその正規化だけから既知 ambient parent membership は出ない。

数学入力は `NIL22_STABLE_PROJECTIVE_ACTION.md` と
`CLASSIFICATION_NIL22_IMAGE_STABLE_CUBIC_20261107.json` の全六磁気源。
新しい CAS、point、operator、rank又は profile を使っていない。
