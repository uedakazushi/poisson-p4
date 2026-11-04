# 二つの squared amplitude の rank2 vector bundle descent

2026-11-04。凍結済み `Q_FIRSTJET_FIVEPOINT_S5_REDUCTION_20261104_v63` の
F=B(s1,s2) に対する別の短い証明。元 packet の bytes は変更しない。
これは独立な2倍率を保持する descent であり、F=B とは言わない。

Y=P(V(2,1))×P(V(1,2)) の一般開とし、L_U,L_V を二つの tautological line bundle
の pullback とする。非零 affine pair の空間は

```
X=L_U^××_Y L_V^×.
```

G0=SL2(A)×SL2(B) の中心 Z=μ2² は、一方の tautological line にだけ−1を
掛ける二つの独立な符号である。したがって square map は厳密に

```
X/Z ≅ {both components nonzero in L_U^⊗2 ⊕ L_V^⊗2}.
```

この rank2 vector bundle は **genuine Aut(Q)=(PGL2×PGL2)⋊C2-linearization** を持つ。
中心 Z は両 square 上で自明であり、追加 Spin5 成分は二つの summand を交換する。
その二乗は G0 の中心なので、square 上では恒等である。よって lift の符号を
選ぶ必要なく、作用は Gext/Z=Aut(Q) を経由する。

五点の一般 configuration に automorphism がないことは元 packet で確認した。
従って対応する Y の Aut(Q)-作用は一般自由であり、適当な quotient open の上で
Y→Y/Aut(Q) は torsor となる。上の genuinely linearized vector bundle はその
torsor に沿って rank2 vector bundle へ descent する。非分裂 ruling や spin lift に
関する gerbe は、square の線形化を作る前提として既に Z を除いているため残らない。

vector bundle の generic fiber は係数体 B=C(Y/Aut(Q)) 上の普通の2次元 vector space。
成分の両方が非零という開条件は関数体を変えないので

```
C(X)^Gext = B(s1,s2).
```

ここで B=C(P(V5))^S5 は元 packet の [3,2] 表現の射影商の関数体である。
これは split な二つの ruling を初めから仮定する証明ではない。二つの line が
交換される場合にも、その**直和 vector bundle** が genuine に descent することを使う。

元 packet の projective no-name により B に4変数を加えると純超越体になる。
従って F=B(s1,s2) に2変数を加えると純超越体となり、stable level は高々2。
B や F 自身の純有理性、全多項式不変式環、NEW25 の式は結論しない。

追加の行列検算、点、jet、GB、有限群不変式探索は行っていない。
