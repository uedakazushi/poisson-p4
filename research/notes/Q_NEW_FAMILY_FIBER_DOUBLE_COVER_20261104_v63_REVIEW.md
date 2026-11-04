# Actual fiberからJ直線への二重被覆：独立監査

2026-11-04。明示actual生成元wの独立レビューを入力とし、QQ/C上のgeneric関数体について
**K_B(w)=K_B(J)(sqrt(r0*J*(J-Jcrit)))** を受理する。
新v0座標(xi0,r0,s0,q0,B0)は旧t=2 sectionの同名座標と区別する。

accepted v0 sectionから `u²=1600r0,K=(8s0+u*q0)/(8s0-u*q0),Iv=u/[512(8s0+u*q0)]`。
前reviewのcompact J(w)へ代入し、sparse numeratorのuべきを
`u^(2j+epsilon)→(1600r0)^j*u^epsilon` と独立にreduceする。
分子差が零、分母が非零であり、保存された新式に一致する。

```
A=200r0*q0*(s0-1), B=s0-25r0*q0²,
D=r0*q0²*(s0²-25r0*q0²), U=1600r0,
J=(256/25)(A+B*w)²/[D(U-w²)]。
```

J D(U-w²)-(256/25)(A+B*w)²の三つの二次係数から、独立に判別式を計算すると

```
disc=(80D)²*r0*J*(J-Jcrit),
Jcrit=256-256/(25r0*q0²)。
```

Jはbase上超越的で、genericにJcrit≠0。
J=0とJ=Jcritでの奇数valuationからr0*J*(J-Jcrit)はbase関数体でも平方でなく、degree2。
二次方程式の完成平方から上のactual field equalityが従う。
square factorを除く際にr0まで平方として落としていない。

分岐点w=-A/Bとw=-UB/Aはgeneric base上有理で、それぞれJ=0,Jcritへ写る。
degree2のrational mapの微分分子は定数倍の(A+B*w)(BU+A*w)なので、これが二つの全分岐点。
極の特殊化やbranchの合流を含む全special fiberを分類したとは言わない。

v0 sectionではR=rf²→256²から `w=u(K-1)/(K+1)=200r0*q0/s0`。
これを上式へ入れるとJ=256を独立に確認する。
この一つの境界sectionからJ=256上の全点や他成分への所属は推論しない。
NEWの定義式・spectral次数・全Poisson成分分類への主張も無い。

両reviewは一つの有限rational-arithmetic checkerで同時に再現できる。

```
python code/q_new_family_explicit_fiber_generator_review_20261104_v63.py \
  --output /tmp/explicit-fiber-review \
  --double-output /tmp/fiber-double-cover-review
```
