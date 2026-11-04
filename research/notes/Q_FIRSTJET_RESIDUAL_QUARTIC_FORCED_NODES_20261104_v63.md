# 平面四次の滑らかさprobe失敗と、m=0上の二つの普遍的節点

2026-11-04。以前のreciprocal quartic packetは変更しない。許可された一つの
**係数特殊化** `(x,y,a,k1,k2)=(1,2,3,5,7)`、素数32003だけを検査した。
beta/source点ではなく、普遍四次の係数点である。Macaulay2の三つのprojective
Jacobian patchesは `r=1:false, m=1:true, z0=1:false`。
infinityの二つのbinary-quartic patchesは両方true。約1.824秒、child RSS93,704KiB。
四次全体のsmoothnessは **失敗として保存**し、別の係数へ移らない。

この結果をgeneric singularityの数値的証拠として昇格するのでもない。
以下は既存の普遍恒等式からの別の厳密な証明である。

## 1. 二つの節点中心は代数的に必然である

reciprocal packetのr,t,m,u,v,c,f,B0,K1を使い、

```
qR=r²+4fuv, qT=t²+4cuv,
D=B0+2ar+2dt,
W=2uvD+(u+v)rt,
4uv*Psi=W²−m²qR qT,      t=(K1−xr)/y.
```

`m=0` でのu=vをu0、他の係数の制限を添字0で表す。
u0はgeneric baseの非零元で、K10はQQ(x,y,a,k2)上超越的である。すると

```
P0(r)=r*t0+u0(B00+2ar+2d*t0),
t0=(K10−xr)/y,
Psi(r,0)=P0(r)².
```

P0のr²係数は−x/y、r係数はK10/y+2u0a−2u0dx/y、定数項は
u0B00+2u0dK10/y。従ってrについてのdiscriminantのK10²係数は1/y²で非零。
P0はgeneric base上で相異なる二根を持つ。この有限etale二重点中心を取り、
必要ならその一つの根alphaを付け、他方をbetaと書く。

各根でqR0が零だとする。rはQQ(x,y,a,k2)上代数的であり、P0=0のK10係数から
r+2u0d=0が必要になる。しかし

```
r²+4f0u0² = 4u0²(d²+f0),
d²+f0 = 2ax+3x²+y² !=0.
```

矛盾する。同様、qT0=0ならt0は同じ係数体上代数的、r=(K10−yt0)/xは超越的。
P0をrの一次式として見るとt0+2u0a=0が必要だが

```
t0²+4c0u0² = 4u0²(a²+c0),
a²+c0 = −2ax−x²+y² !=0
```

となり矛盾する。従って **両中心でqR0*qT0は単元**。

またP0のr微分が各根で非零なので、局所座標(m,W)が取れる。
式はunit倍の `W²−m²*(unit)`。標数0の代数閉体上で二次初項は非退化なので、
これは二つの相異なる滑らかな枝を持つ **通常の節点**である。
よって普遍平面四次は、generic base上でm=0に正確に指定された二つの通常節点を持つ。
「他には特異点が無い」ことはここでは証明しない。

今回の同じ係数点では

```
P0(r)=−(2r²−287r+10914)/4,
disc(2r²−287r+10914)=−4943 !=0 mod32003.
```

これが元のsmoothness testの失敗を構造的に説明する。数値的失敗を理由に
universal quartic identityや粗い種数上界17を撤回する必要はない。

## 2. 各幾何的四次成分の正規化種数は高々1

平面四次が幾何的に既約なら、算術種数3から上の二節点のdelta各1を引き、
正規化種数は高々1。可約またはnonreducedなら、各幾何的既約support成分の次数は
高々3なので、標準plane genus boundからやはり高々1である。
従って **全ての幾何的成分について商曲線の正規化種数は高々1**。
滑らかな楕円曲線であること、全四次の既約性、full二次coverの正確な種数は未決である。

## 3. 別模型への具体的な次の操作

alpha,betaを二節点のr座標とし

```
U=r−alpha*z0, V=r−beta*z0, M=m,
delta=alpha−beta,
r=(alpha*V−beta*U)/delta,
z0=(V−U)/delta.
```

元の四次の両点での消失次数が2なので、変換後のquarticはU,Vの次数がそれぞれ
高々2になる。M=1へdehomogenizeした式をH(U,V)とすれば、これは(2,2)の式である。
P¹_U×P¹_Vへの閉包は、二節点をblowupして両者を結ぶ直線をcontractした模型に対応する。
四次のstrict transformとその直線のintersectionは4−2−2=0なので、このcontractは
曲線の近傍では同型である。

次に許可されたdistinct methodは、この **(2,2) proper model** を同じ係数点で一度だけ
検査することである。新しいbeta/source点、別の係数、元四次smoothnessのretryではない。
alphaは上の二次etale algebraで保持でき、二根を別々に数値選択する必要もない。

全四patchでJacobianが空なら、proper flat (2,2) familyとsmooth locusのopen性から
generic smoothnessへ進める。smooth (2,2) divisorは幾何的にconnectedである
（非交叉のeffective成分は両方向を合わせて(2,2)にできない）ため、geometric
irreducibilityとgenus1も従う。ただし実際の新jobの結果は別receiptで記録する。

## 保存

失敗した一回のsource:
`code/q_firstjet_residual_quartic_smooth_specialization_20261104_v63.py`。
全15特殊係数、M2 script/log、false report:
`reports/Q_FIRSTJET_RESIDUAL_QUARTIC_SMOOTH_SPECIALIZATION_20261104_v63/`。
一回のM2は正常終了したが、数学的なsmoothness gateは通らなかった。
このnoteの二節点証明はgeneric coefficientによる厳密な構造補題であり、別点やGBを使わない。
