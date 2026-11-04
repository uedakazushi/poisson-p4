# Universal ordered residual curve：四つの無限遠分岐、種数3、m次数8

2026-11-04。対象はfull W16 imageのgeneric spectral fiberを、既存のordered-eigenvalue
baseと相対符号quotientの後で書いた曲線である。NEW第一jet hypersurfaceの有限
multisection自体や、その定義式ではない。以下は新しい係数点・beta/source点・不変式
jobを用いず、proper (2,2)模型の受理を前提に分岐を全て数える証明である。
一次報告はcombined independent review待ちとして保存する。

## 1. 既知の模型を固定する

```
K=QQ(x,y,a,k1,k2),  B=QQ(x,y,a,k2),
qR=r²+4fuv, qT=t²+4cuv,
t=(K1−xr)/y, D=B0+2ar+2dt,
W=2uvD+(u+v)rt,
4uv*Psi=W²−m²qR*qT.
```

係数c,f,u,vはm一次、B0はm二次、K1=k1/2−A0(m)である。全実係数はreciprocal
packetで固定し、P0 quadratic node packetとproper (2,2)模型は同じ係数を使う。
そのproper模型がgeneric geometrically smooth・connected genus1を与えるため、
平面四次の正規化Gammaは幾何的に既約な種数1曲線である。
元の平面四次の特異点はm=0上の受理済み二節点だけ。

full curveの関数体候補は `K(Gamma)(h), h²=qR`。K(m)上の全generic intersection
algebraは既受理の四次×二次表示により次元8である。

## 2. 有限部で qR=qT=W=0 は存在しない

この三式はk1を含まない。その共通零点が存在すればmは以下の非零norm polynomialの
根になる。

`q=uv, l=u+v` と書き、free rank4 algebra

```
B[m][r,t]/(r²+4fq, t²+4cq)
```

でWのmultiplication determinantを取る。二回のquadratic norm、または明示的4×4行列から

```
E=B0²−16fq*a²+16cq*d²−4cf*l²,
F=a*B0+2cd*l,
Norm(W)=16q⁴*(E²+64fq*F²).                     (1)
```

専用の小さい普遍QQ演算で4×4 determinantの全係数と(1)を比較した。
実際のm-leading coefficientsは

```
u1=−(x−y)²/(4xy), v1=−(x+y)²/(4xy),
q2=Delta²/(16x²y²), l1=−(x²+y²)/(2xy),
B0,2=2, c1=f1=1.
```

従ってEのm四次係数は

```
4−4l1²=−Delta²/(x²y²),
```

第二項64fqF²はm次数高々7なので、(1)の全m次数とleading coefficientは

```
degree_m Norm(W)=16,
LC_m Norm(W)=Delta^12/(4096*x^12*y^12) !=0.      (2)
```

これを巨大な全base係数に展開する必要はない。Normはfree monic quadratic algebra二段の
determinantなので、共通零点は必ずその零点へ送られる。従ってmはB上代数的、さらにr,tも
それぞれmonic quadraticの根としてB上代数的になる。
ところが `xr+yt=K1=k1/2−A0(m)` は超越元k1をB上代数的にしてしまう。矛盾。
よってgeneric Kの代数閉体上で **有限 qR=qT=W=0 は無い**。

## 3. 分母boundaryを残した有限零点の偶性

Gammaの有限点でqR=0を考える。m=0上では曲線の点は二つの節点の枝だけであり、
両節点でqRは単元であることを既証。正規化上の二枝でも単元のままなので、ここに零点は無い。

u=0またはv=0ならmはB上代数的で、genericにm*fと他方のu/vは非零。
qR=r²だからr=0。式(3)のPsiを直接使えば

```
Psi|_(uv=0,r=0)=−m²*f*t²,
```

従ってt=0となり、K1(m)=0がk1の超越性に反する。u/vの同時零やfとの同時零は
k2の独立性に反する特殊base条件であり、generic Kではない。

参考としてf=0も直接検査できる。m,u,vはgenericに非零、qR=r²からr=0となり、
`Psi=uv*(B0+2dt)²`。dはBの非零元、tはk1について非定数なので零にはならない。
c=0は分岐判定に分母として使わない。c=0を含め、次の恒等式はmuv以外の反転を要求しない。

残るqR零点ではmuvは単元。恒等式からW=0であり、§2よりqTは非零、すなわち局所単元。
従ってGammaの正則局所環で

```
qR=(W/m)² / qT.
```

**全ての有限零点のvaluationは偶数**である。affine polynomialであるqRに有限poleは無い。
有限部には二次coverの分岐点は存在しない。

## 4. 無限遠の四点では pole order が正確に3

一回の旧係数probeに含まれていたbinary quartic at infinityの二patch検査は双方PASSだった。
この部分のproper familyはP1上のdegree4 divisorで、r⁴係数16x⁴y⁸が単元。
finite flat degree4のfamilyでreduced locusはopenなので、generic infinity divisorは
**相異なる四つの幾何点**である。これは元の平面四次全体がsmoothという失敗したgateとは別。

無限遠ではm=0にできない。実際、homogeneous quarticの `[r:0:0]` での値は
16x⁴y⁸*r⁴で非零。従って各無限遠点でmは単純pole、r/mは有限である。
qR=r²+4fuvの全次数3 leading partは

```
Delta²/(4x²y²) * m³.
```

係数は非零で、r²項はpole order高々2なので相殺できない。よって四つの点すべてで

```
ord(qR)=−3.
```

二つの節点は有限部なので、proper normalization/blowupは無限遠を変更しない。
標数0の平方根coverは、この四つの奇数poleでそれぞれ単純分岐する。
§3と合わせて **分岐divisorの次数は正確に4**。

## 5. Connectedness、種数3、actual selected componentの次数8

奇数poleを持つqRは幾何的関数体でも平方ではない。従って二次coverは
geometrically connected・irreducibleで、滑らかなproper normalizationをCとすると

```
2g(C)−2=2*(2g(Gamma)−2)+4,
g(Gamma)=1 なので g(C)=3.
```

これは幾何的種数であり、Prymとの同定、hyperellipticity、元のWeyl quotientの種数は主張しない。

四次商がgeometrically irreducible、さらにこの二次coverもgeometrically irreducibleなので、
元の全generic intersection algebraはK(m)上の **field** でdimension8。
既受理の実際のW16像成分はmが非定数で、このalgebraのfield factorを与える。
fieldに非自明なfield factorは無いから、そのselected m-dominant成分は全体を占める。
従ってこの明示したordered-base/relative-sign chartの後で、

```
actual selected componentのm-degreeは正確に8、正規化種数は3。
```

NEW first-jet hypersurfaceはこのuniversal曲線上の別のfinite multisectionであり、ここでの8を
そのspectral degreeへ読み替えない。既知のNEW target 10/20/40とは対象が異なる。
元のbaseでのWeyl descent、全boundary fiber、NEW定義式、Poisson component分類は未決のまま。

## 保存

専用source `code/q_firstjet_residual_genus3_branch_20261104_v63.py` は、4×4 universal normの
全係数照合とgate/hash確認だけを行う。大きいGB、別係数点、beta/source点、不変式jobは無い。
`reports/Q_FIRSTJET_RESIDUAL_GENUS3_BRANCH_20261104_v63/` にmatrix、branch schema、scopeを保存。
proper (2,2)模型と全branch論証のcombined independent reviewを受けてから台帳へ採用する。
