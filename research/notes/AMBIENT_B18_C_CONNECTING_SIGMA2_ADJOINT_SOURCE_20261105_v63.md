# 残存 sigma2 の実三行 adjoint と四次 primitive frontier

残った C*sigma2_rr の係数行 C は列0,16,20にのみ非零である。
既存の固定24単元逆から n0 を再構成し、L=C*n0 を取ると、
行8,12,13にのみ係数 −43303/59400, −149887/237600, 30713/59400 が残る。
これは B sector の z0,z4,z5 Hadamard の三行であり、r³で割る行を含まない。

固定 section n=Esel*(mEsel)^{-1} に対する cycle の微分は
sigma1=−n0*m1*sigma0、sigma2=−n0*(m2*sigma0+m1*sigma1) である。
n の微分を落とした式ではない。したがって C*sigma2 はこの三行の m2 と
m1*sigma1 の contraction だけに依存する。

ordered Hadamard を元の z=0 上で評価する。適応座標 y_i=z_i+c_i r² の原点
y=0 と取り違えない。必要な m2 の次数は実 B 関数の三次までであり、
最高次は z0*r²,z4*r²,z5*r²。既存の beta2 と二次係数による prefix carry も残す。
実 sigma0 で全56列を contraction すると、五つの B 関数のうち一つは零、
一つは nuA*bA67 で既存低次項のみ、残る三つは delta18/60/67,delta63,bB63 の
固定線形結合になる。

別途、V4*sigma0 自体の未知四次項も全56列の source から列挙する。
first7 bAi=zAi と、bA63|z=0=gA の実 telescope を用いると、未知関数は
delta18,delta60,delta67,delta72,delta63,gA,bB63 の七つに収まる。
各々の r³ell,r²ell²、計十四係数の上界である。tauAB からの既知の
−nuB2*bB2 の四次 carry は別に保持する。

これは V4*sigma0 部分だけの frontier であり、接続写像全体の最高次依存ではない。
V3*sigma1、混合 V2*sigma1、L*m2、L*m1*sigma1 を合わせた coherent contraction
が次の SOURCE gate である。sigma1_r の m1 最終行には既存四次 source が
関わり得るため、先に零とはしない。既存 joint105/選択135/九第四 adjoint や
二族の一変数六次 cache は、全 sigma2 や全平面四次表とは異なる。
本証明は新しい高次係数や接続写像の値を一つも問い合わせない。
