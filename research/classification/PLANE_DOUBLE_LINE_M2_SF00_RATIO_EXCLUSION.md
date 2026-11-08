# R=2L+M・m二SF00：唯一模型の全higher源はratio-Casimirで空

2026-11-08、root純手証明、独立確認前。計算実行無し。
`PLANE_DOUBLE_LINE_M2_SQUAREFREE_WHOLE_FIRST_SYMBOL_REDUCTION.md` の全必要還元を使う。
その00唯一模型はf=a²b、X=H−2a∂c、Y=−H−2a∂c、H=diag(1,-2,1)、
normal K=kappa*u*v、kappa=−5、actualG2=−27a²buvである。
原稿で残していたこのfullhigher義務も、以下の最弱構造証明で閉じる。

全raw36 P1と全P2九を許す。row divergence Ru=kappa*u、Rv=−kappa*v。
canonical divergencezero strict liftPi=P+E∧R/5はactualFをCasimirに持つ。
zeta=v/u、K=C(u,v)としてM=Pi sharp d(zeta)を取ると、全sourceで

```
M(zeta)=0,
M(u,v)=(3*k_eff/5)*(u,v), k_eff=kappa*zeta=−5*zeta,
M(base)=(A(zeta)−2*k_eff*I/5)*base+B(u,v),
A(zeta)=−zeta*(X−Y)=−2*zeta*H.
```

Bは全higher mixed coefficientのrational normal-homogeneous degree一vector。
base-base P1/P2はd(zeta)contractionへ寄与しない。
Mはdense normal fieldで読む恒等式であり、u/v又はunknown source coefficientを全sourceの仮定として反転しない。

Aの固有値はs,-2s,s、s=−2zeta。k_eff/s=5/2である。
translation equationのdetは(k_eff-s)²(k_eff+2s)で、Kのunit。
T=−(k_eff I−A)^−1B、base'=base+Tによって全translationが消える。
逆base=base'−Tがあり、normal total degreeは一のまま。
この有理translationをglobal PGLとは呼ばない。

全degree五Fをbase degree d、normal total degree m=5−dで分ける。
元normal order二なのでm0/1は零で、translation後も生じない。
m3/4/5の全operatorはSym^(5−m)(A)+k_eff(m−2)I。
base spectrum s*(d−3j)を全て入れると、その非可逆比の和は

`k_eff/s∈{-2,1,4,-1/2,0}`。

translation失敗もこの和に入る。5/2は含まれないので、全m3/4/5 summandsは零。
従って全Fは元leadingそのものを有理base translationで戻した

`F=−27*u*v*(a+T_a)²*(b+T_b)`

に等しい。K[a,b,c]で正のbase degreeを持つ因子を二つ以上持つ。
元FはC[u,v,a,b,c]上polynomialであるから、Gauss lemmaにより既約性に反する。
全P2九と任意のinhomogeneous particularを保った全source空性である。

この証明と全first-symbol還元及び01/02/03/11/12全空を合わせれば、
**R=a²b・m二generic SFの全被約既約quintic sourceは空**。
そのwhole幾何必要還元と本ratio証明の独立受理前に台帳のG/Eを完了としない。
normaldouble/m三/全分類は別義務である。
