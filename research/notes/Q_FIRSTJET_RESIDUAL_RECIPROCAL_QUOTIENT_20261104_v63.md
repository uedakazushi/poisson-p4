# 固定 spectral 残余曲線の reciprocal involution と平面四次商

2026-11-04。既存の二つの (2,2) Laurent 方程式から、全 generic intersection を
**平面四次曲線の二次 cover**へ有理的に書き換える。新しい beta 点、source 点、
不変式展開、GB、八次式の因数探索は用いない。

ここで得るのは ordered-eigenvalue base と既存の相対符号 quotient の後の模型である。
各幾何的 m-dominant 成分の正規化種数は粗い上界17を持つが、正確な種数、四次商の
滑らかさ・既約性、選択された成分の m-degree8 はこのpacketでは証明しない。

## 1. 元の二式と generic base

凍結済み `Q_FIRSTJET_FIXED_BASE_RESIDUAL_CURVE` と
`Q_FIRSTJET_RESIDUAL_OCTIC_PRESENTATION` の記号を使う。

```
K=QQ(x,y,a,k1,k2),             mはK上超越的,
Delta=x²−y², S=x²+y², ell=(a+x)²−y²,
d=−(ax+S)/y, c=m−ell, f=m−x²ell/y²,
B2=−4a²x²−2ax³−6axy²+3x⁴+y⁴,
u=(k2−B2−(x−y)²m)/(4xy),
v=(k2−B2−(x+y)²m)/(4xy),       u−v=m.
```

さらに

```
A0=ad(ay+dx)+(ay+dx)(u+v)+cdy+fax,
K1=k1/2−A0,
B0=(a²+c)(d²+f)+m²+2ad(u+v)−Delta².
```

generic fieldでは x,y,Delta,m,c,f,u,v は非零。K1は
`QQ(x,y,a,k2,m)` 上の超越元である。対象は

```
L1=−cyT−xvR+fxu/R+yuv/T−K1,
L2=B0+cRT+fu²/(RT)−cfT/R−v²R/T
           −2cdT+2afu/R−2avR+2duv/T.
```

元の全 generic scheme algebraは既受理の八次表示により K(m) 上dimension8。
special mや表示分母の消えるboundaryを消去したglobal theoremとは区別する。

## 2. 同時 reciprocal involution

次はmを固定する有理involutionであり、**両式を同時に保つ**。

```
iota(R,T)=(-fu/(vR), -uv/(cT)).                  (1)
```

商の座標を新しい小文字 r,t,z で

```
r=−vR+fu/R,
t=−cT+uv/T,
z=(vR+fu/R)(cT+uv/T)
```

と置く。r,t,zは(1)で固定され、

```
z²=(r²+4fuv)(t²+4cuv),
L1=xr+yt−K1,
L2=B0+2ar+2dt+((u+v)rt+mz)/(2uv).               (2)
```

全て普遍的なLaurent polynomial identityであり、標本の一致ではない。

このinvolutionの行列による説明もある。元のadjoint対を同時に
`(C,D)→(−C,−D)` とし、`G=diag(J2,J2)`、`J2=[[0,1],[-1,0]]` で共役すると
Dの順序付き対角形が戻る。次に `diag(s,s^-1,t0,t0^-1)`、
`s²=1/c,t0²=1/f` でb=e=1を戻せば、off-diagonal係数は

```
g' =−(s/t0)j, h'=s*t0*i,
i' =h/(s*t0), j'=−(t0/s)g.
```

従ってR'=g'/h',T'=g'h'は正確に(1)。s,t0を選ぶ有限拡大があっても、
誘導されたR,T式には平方根が残らない。actual quadratic covariants C,D の同時負号は
betaのi倍に対応するため、元の斉次像全体を保つ。ただしordered finite base-change後の
個々の選択成分を保つかどうかを、この事実だけで仮定しない。

## 3. 四次商の明示式と実際の degree cancellation

(2)から

```
t=(K1−xr)/y,
D0=B0+2ar+2dt,
z=−(2uvD0+(u+v)rt)/m.
```

次を置く。

```
Psi = r²t²+(u+v)rtD0+uvD0²
                      −m²(c r²+f t²+4cfuv).     (3)
```

恒等式

```
(2uvD0+(u+v)rt)²
 −m²(r²+4fuv)(t²+4cuv) = 4uv*Psi               (4)
```

により商の式はPsi=0となる。rについては明らかに四次で、最高係数はx²/y²。
重要なのは、§1の実際のc,d,f,u,v,K1,B0を代入すると、見かけのm六次以上が消え、
**(r,m)の全次数が正確に4**になることである。

一回の小さいQQ全係数計算でこれを検証した。名前付き分母 `(4xy)^4 y^8` を払い、
共通monomial x²y²とcontent16だけを除いた整数係数四次式は

```
F(r,m)=16x²y^10 * Psi(r,m).
```

Fの全15個のr^i m^j係数（i+j<=4）は全て非零で、係数変数(x,y,a,k1,k2)を
含めた項数は1592。全係数表を保存する。四次のr⁴係数は16x⁴y⁸。
従って F はP²_(r:m:z0)の平面四次を与える。ここでz0はhomogenizing coordinateで、
§2のzとは別物である。

## 4. 全八次 algebra の quartic-times-quadratic 表示

`q=r²+4fuv`、`h=vR+fu/R` とする。次の式が全generic schemeの有理逆を与える。

```
h²=q,
R=(h−r)/(2v),
T=(z/h−t)/(2c).                                (5)
```

必要なhの単元性もgeneric fieldで証明できる。Psi=q=0の共通根を仮定する。
rは `K00=QQ(x,y,a,k2,m)` 上代数的である一方、t=(K1−xr)/yはK00(r)上超越的。
(4)より `2uvD0+(u+v)rt=0`。これはtの一次式なので、その定数係数から
`B0+2ar=0` が必要となる。q=0と合わせると

```
B0²+16a²fuv=0.
```

しかし左辺はm四次で最高係数4。全実係数を代入してこの非零最高係数も検証した。
よってqは四次algebraの単元で、hも単元となる。

また `(h−r)(h+r)=4fuv`、
`(z/h−t)(z/h+t)=4cuv` は単元なので、(5)のR,Tも単元である。
これらの式と(2)–(4)を使えば、両方向の式は元のr,t,zおよびR,Tを回収する。
従って **零点だけでなく全algebra** について

```
K(m)[R^±1,T^±1]/(L1,L2)
  ≅ K(m)[r,h]/(Psi, h²−q).                     (6)
```

右辺は四次商の上にmonic二次式を重ねたfree rank2 algebraであり、全次元8。
iotaはh→−h、その固定algebraは正確にK(m)[r]/(Psi)。qが単元なのでこのgeneric
二次coverはetaleであるが、mの特殊値やprojective boundaryでの分岐が無いとは言わない。
既約性や被約性を仮定しない全scheme algebraの表示である。

## 5. 無条件の粗い種数上界17

係数体Kを代数閉体へ拡大し、mに優越する任意の幾何的成分を取る。
(6)により、その曲線の関数体は、平面四次Fのある既約成分の関数体から
qの平方根を加えたもの、または二次式がsplitした場合の同じ関数体である。

平面成分の次数をd<=4、正規化をGammaとすると、標数0のplane genus boundから

```
g(Gamma)<=(d−1)(d−2)/2<=3.
```

q=r²+4fuvはaffine(r,m)の全次数3以下の多項式。Gamma上のpole divisorの次数は
高々3d<=12である。これはhomogeneous cubic numeratorとline-at-infinity³の
divisorを比較すればよく、境界での相殺はpoleを減らすだけである。
零divisorの次数も同じなので、奇数valuationを持つ点の総数は高々24。
非split二次coverの正規化をCとすると、標数0のRiemann–Hurwitzにより

```
g(C)=2g(Gamma)−1+(ramification degree)/2<=17.     (7)
```

splitならg(C)<=3。これで **ordered-eigenvalue拡大と既存相対符号quotientの後の、
任意のm-dominant幾何成分の正規化種数は高々17** と証明される。
既受理の実像成分はmが非定数なので、この範囲に入る。
これは正確な種数17、平面四次の滑らかさ、四次商全体が一成分、全八点が同一成分、
元のWeyl quotient前後で同じ種数という主張ではない。

## 6. 正確な種数のための次の有限gate

次に調べる対象は元の三変数complete intersectionのnaive genus33ではなく、
保存済みの一つの平面四次と一つの二次coverである。

1. Fのhomogenizationと三つの偏微分のprojective同時零点を調べ、generic smoothness
   とgeometric irreducibilityを証明する。係数のある一つの良いspecializationを使う場合は、
   全係数がp-integralなproper平面四次familyとsmooth locusのopen性を明示する必要がある。
2. 正規化上でqの全零・pole divisorの奇偶を求める。finite boundary
   `m*c*f*u*v=0` とinfinity、qと `t²+4cuv` の同時零も保持する。
3. 二次coverがconnectedか、選択されたactual componentがどの因子に当たるか、
   finite ordered-base/Weyl quotientの識別を別々に確認する。

例えば四次が滑らかで、分岐点がinfinityの四点だけならRiemann–Hurwitzは種数7を与える。
しかしその仮定は本packetでは証明しておらず、**種数7を結果として採用しない**。
任意の高次数invariantやNEW first-jet hypersurfaceの定義式・次数をここから推測しない。

## 7. 保存と実行

専用sourceは `code/q_firstjet_residual_reciprocal_quotient_20261104_v63.py`。
同名reportsに普遍identity、全15係数、1592項、runtimeと入力hashを保存した。
実測0.221472秒、最大RSS51,648KiB。60秒/2GBの一回のgeneric symbolic operationであり、
新しいbeta/source点・finite-field probe・GB・octic factorizationはない。

```
.venv/bin/python code/q_firstjet_residual_reciprocal_quotient_20261104_v63.py \
  --output /tmp/reciprocal-quotient-replay
```

上界17は保存したdegree4/degree3と標準curve論からの別の証明であり、正確な種数を
計算したCAS出力として報告しない。別のGENUS_BOUND receiptがこの依存と範囲を記録する。
