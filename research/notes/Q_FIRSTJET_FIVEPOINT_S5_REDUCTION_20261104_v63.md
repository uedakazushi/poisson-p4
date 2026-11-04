# 五交点・S5 の5次元表現への還元と残る二倍率

2026-11-04。対象は複素数体上の W16 の Sp4 不変式体 F。
既受理の Spin4 biform 切片と追加 normal-line 成分を用いる。
新しい finite-group invariant search、β 点、高次数評価、GB は行わない。
**P4/S5 や F の純有理性は主張しない。**

## 1. 五点から各形式を一意に戻す

固定した滑らかな Q=P1×P1 上の U∈H0(O(2,1)), V∈H0(O(1,2)) を、
それぞれの非零定数倍を除いて考える。一般の二形式は横断的に交わり、
Z=U∩V は五つの相異なる点となる。Koszul 列を O(2,1) で twist すると

```
0→O(-1,-2)→O ⊕ O(1,-1)→I_Z(2,1)→0.
```

Künneth により H0(O(1,-1))=H0(O(-1,-2))=H1(O(-1,-2))=0。
従って h0(I_Z(2,1))=1。同様に h0(I_Z(1,2))=1。
よって Z は各形式を定数倍を除いて一意に決める。

P5×P5 と五点の一般 configuration 空間は共に既約10次元。
交点写像は上の開集合で一対一なので、この写像は五点 configuration の稠密開へ
双有理である。ruling を交換すれば二つの形式も交換され、追加成分の符号は
ここでは projectivization により消える。従って projective pair の全商は
五つの unordered points の Aut(Q)=(PGL2×PGL2)⋊C2 商に双有理。

## 2. P3 の projective frame と quadric の5次元表現

|O_Q(1,1)| で Q を P3 に埋め込む。五つの一般の ordered points は projective frame。
その frame を、V4 の五ベクトル f0,…,f4（f0,…,f3 は標準基底、f4=−Σf_i）の
五つの projective lines に一意に移せる。

この五点を通る quadratic forms は10係数に独立な五条件を課し、5次元線形空間 V5
をなす。一般の quadric は非退化。例えば x0x1−x2x3 は全五点を通り、対称行列の
determinant は1/16。この例は quadric 線形系の非空性だけで、β の点評価ではない。

frame の relabeling は S5 であり、f_i を置換する標準4次元表現に正確に lift する。
従って、configuration と quadric の quotient をどちらから取っても

```
B := C(projective pair)^Aut(Q) ≅ C(P(V5))^S5.
```

V5 は Sym²(Std4*)→permutation5（q↦(q(f_i))) の核。
基底を六つの off-diagonal monomial の最初五つから最後 x2x3 を引いたものとする。
隣接 transposition 四つの5×5行列を保存し、全 Coxeter 関係と全120像の
projective faithfulness を確認する。conjugacy class 順
(1⁵),(2,1³),(2²,1),(3,1²),(3,2),(4,1),(5) の character は
(5,1,1,-1,1,-1,0)、内積1。従って V5 は5次元既約 [3,2] 表現である。
ここで有限群不変式の探索はしていない。

## 3. ruling orientation と二つの scale を失わない

F は projective pair の商 B ではない。U,V の二倍率を戻す必要がある。
H=Gm_U×Gm_V とし、連結群 SL2(A)×SL2(B) の中心 (ε_A,ε_B) を見ると、
U,V の倍率は (ε_B,ε_A)。従って四つの中心はちょうど独立な二つの符号。
effective amplitude torus は

```
H'=H/μ2²,  X*(H')=2Z e_U ⊕ 2Z e_V.
```

つまり U²,V² の倍率を使う。これらは PGL2×PGL2 の genuine な表現に descend する。
追加成分は二つの effective character を交換する。

B 上の一般 quadric の二 ruling は、必ずしも個々に B 上で定義されない。
その quadratic étale algebra を E/B と書く（split の可能性も許す）。
frame chart では cover は ε²=det(q) で、q を引き戻す standard4 行列 g に対し
ε は det(g) 倍される。S5 ではその符号を保持する。q の projective scale λ に
対して ε は λ²倍となる。従って ruling label を単なる無条件の S5-sign だけや、
初めから split な二成分として扱わない。

## 4. squared amplitude が gerbe の障害を消す理由

一般 configuration は Aut(Q) stabilizer を持たない。これは V5 の projective
S5 作用が faithful であり、非恒等元の固定 locus が有限個の真の閉集合だからである。
したがって B の separable closure 上で一つの affine pair を選び、Galois transport と
比較すると、projective pair を結ぶ Aut(Q) の元は一意である。

それを SL2×SL2（および ruling 交換成分）へ lift した cocycle の失敗は、
中心 μ2² だけに値を持つ。amplitude を H' へ写すとこの失敗は零になる。
残る倍率は、二 ruling を交換する monodromy によって twist された torus の
**真の1-cocycle**である。その torus の character lattice は上の二つを単に交換するので

```
T = Res_(E/B) Gm.
```

従って F の B 上の一般 fiber は T-torsor である。非分裂 quadric の ruling conic が
P1 でなくても、U,V の square は中心 character を殺しているので、追加の spin-center
または Brauer の2-cocycle は残らない。これは unsquared amplitudes のままで
torsor が自明と仮定する議論ではない。

Hilbert90（E 上の1次元ベクトル空間が自由であること）により H1(B,T)=0。
よってこの torsor は自明。T は2次元ベクトル空間 Res_(E/B) A1 の norm 非零開部分なので
B 上有理であり、非標準的に

```
F ≅ B(s1,s2),  B=C(P(V5))^S5.
```

を得る。二つの independent scale を数え落とした F=B ではない。

## 5. projective no-name と stable rationality の改善

W4=Std4 とする。S5 は P(V5) と W4 の双方で一般自由に作用する。
Y=(P(V5)×W4)/S5 を二通りに見る。

B=P(V5)/S5 の上では Y は associated rank4 vector bundle なので、

```
C(Y)=B(t1,t2,t3,t4).
```

一方 W4/S5 は rational4次元で、五変数の sum=0 を使えば、その関数体は
基本対称式 C(e2,e3,e4,e5)。この上で Y は**線形5次元表現 V5 に付随する
rank5 vector bundle の projectivization**である。従って projective generic fiber は
Brauer-twisted ではない P4、

```
C(Y)=C(e2,e3,e4,e5)(u1,u2,u3,u4)
```

は8変数の有理関数体となる。これが projective no-name による B の stable level 上界4。
affine V5 を使って余分な一変数を数え込んだ議論ではない。

F=B(s1,s2) だから、F にさらに二変数を加えると B に四変数を加えた体となる。
従って **F の stable rationality level は高々2**。
F または B 自身の純有理性、stable level の最小値は依然主張しない。

## 6. 計算と適用範囲

```
PYTHONDONTWRITEBYTECODE=1 .venv/bin/python code/q_firstjet_fivepoint_S5_reduction_20261104_v63.py
```

許可された一回の小計算は quadratic space、四生成行列、Coxeter、character と
120個の projective faithful 像だけを確認する。生成不変式の探索や新しい β 点は使わない。
幾何的対応、center・ruling の torsor descent、Hilbert90、projective no-name は
上の証明であり、行列の標本から推測したものではない。

NEW25 の未知 divisor、次数24の式、残余曲線 genus、多項式不変式環の全生成は
今回の結論ではない。凍結済みの既存原稿・証明書は変更しない。
