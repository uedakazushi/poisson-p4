# concurrent m二doubleの残二A8源は各一実軌道閉包

2026-11-08、PRIMARY E。前段ONE三全source batchの後の純actual normalizer導出である。新しいM1、点、ambient profile、minor又はjetを計算しない。以下の恒等式の独立literal監査は別義務とする。

## 1. 対象と全基底

`CLASSIFICATION_PLANE_M2_CONCURRENT_DOUBLE_THREE_WHOLE_SUPPORT_20261108.json` のh3,b=−2二源を用いる。period tは5/2又は−1/2、E=−2/(t+2)は−4/9又は−4/3、k=6Eである。元座標は(a,b,c,u,v)。全P0は

```
f=ab(a−b), Y=(ta−2b)∂c, X=Y+(−3a+6b)∂c,
Z=(a+3Eb)∂b−3Ec∂c,
P0=Jac(f)+uX∧∂u+(uZ+vY)∧∂v+6E u²∂u∧∂v.
```

K1の全八free coordinatesは、保存basis順で

```
p0=B_c, p1=C_c, (p2,p3,p4)=D, (p5,p6,p7)=A,
```

つまりP_cu[uv]、P_cu[v²]、全P_iv[u²]、全P_iv[uv]である。他の係数は保存完全kernelの列で全て決まる。FULL9P2は全P1値で一意に解け、元全350式は恒等零なので、両源は実Poisson A8族である。m=p1零ではactualDにu因子がある。

## 2. 実六shearで全cohomologyを二係数にする

base-from-normalの有限shear a_i′=a_i+s_i u+t_i vを用いる。P0へのweight一作用をcanonical Euler gaugeで戻すと、全free coefficientsへの作用は

```
B=−DΔ t,   C=0,   D=L s,   A=DΔ s+L t,
Δ=(−3a+6b)∂c,
L=6E I−DZ=[[6E,0,0],[−1,3E,0],[0,0,9E]].
```

det L=162E³は両源で非零定数である。従って任意の元kernel係数に対し

```
s=L^{-1}D, t=L^{-1}(A−DΔ s),
r=p0−(3t_a−6t_b), m=p1
```

が全域定数可逆なshear/residual座標となる。逆shearでP1=rK0+mK1へ移せる。全九P2へのfinite作用はidentity-linearのtriangular translationであり、全方向を保持する。さらにJ2がinjectiveであるため、実変換後のP2はこの新P1に対して保存した一意particularと一致する。

このnormalized二parameter sourceのP2は、t=5/2ではPab=−16rm v²/243、t=−1/2ではPab=16rm v²/81のみである。他八係数は零。

## 3. P0を保存する実unipotentがrを消す

全εについて次のactual GL5 mapを用いる。

```
a′=a, b′=b, u′=u,
v′=v+εu,
c′=c+ε[L_a a+L_b b],
L_a=(3E+1)/(3E²), L_b=−1/E.
```

逆はεを負にするだけである。このmapはP0を全て保存する。実際、base Jac(f)はc-translationに不変、X/Yはc方向で不変であり、mixed Zの追加項は

```
(DZ_AB+3E I)(L_a a+L_b b)+Δ_c=0
```

によって消える。normal bracketも不変である。両periodではL_a=−9/16、L_bは9/4又は3/4。

K1のbad-vectorによるB_a,B_bをm(A_m,B_m)と書くと、

```
(A_m,B_m)=(16/27,10/27) 又は (−16/9,−2/3),
L_a A_m+L_b B_m=1/2.
```

normalized(r,m)をこのmapで押し、canonical Euler gaugeへ戻した全free coefficientsは

```
p0′=r−(5/2)εm, p1′=m,
D_a′=−ε² A_m m, D_b′=−ε² B_m m,
D_c′=−ε²r+(ε³/2)m,
A_a′=2ε A_m m, A_b′=2ε B_m m,
A_c′=2εr−2ε²m.
```

全P2はgradingを保つmapからP1へ寄与しない。再度§2の六shearでD/Aを消した後、rの補正は

```
3t_a−6t_b=εm[A_m(3E−2)/(3E²)−4B_m/E]=0
```

である。従って**全cohomology上の正確な作用は(r,m)→(r−5εm/2,m)**。

m≠0ではε=2r/(5m)でrを零へ移せる。m=0は元のactualu因子補集合として既に保持しており、この反転で落としていない。全finite mapsとshearは明示逆を持つ。

## 4. mのactual scaleと、各一canonical軌道

diag(1,1,1,λ,λ)はP0を保存し、全P1係数r,mをλ^{-1}倍、P2をλ^{-2}倍する。r=0上、λ=m≠0でmを一にする。唯一P2 particularはr=0で零だから、全m非零の点は

```
Pcanonical=P0+K1（bad-vector列一）
```

という同一実PGL軌道に入る。これは恣意的なanchorの選択ではなく、全A8 sourceへの明示逆normalizerである。両periodのcanonicalは異なる可能性を残し、同一性や既親の名前を仮定しない。

各canonicalの実Dは、保存した全Dをp1=1、他係数零へ戻して得る。cについての二次式であり、c²係数はそれぞれ−256u³/243、256u³/9。定数項のv四係数が非零一次式なので係数全体はprimitiveである。c-discriminantをv=0へ制限すると、それぞれ−512 a³u⁵/81、512 a³u⁵となる。非零の非平方monomialなので全discriminantはC(a,b,u,v)上で非平方。従って両Dは絶対既約であり、実rank四点が存在する。

m非零はA8にdenseであり、その全像が一軌道なので、全PGL×A8像閉包はこのcanonical軌道閉包に等しい。m=0の全sourceもこの閉包に入る。従って後続のcanonicalの既知親／全成分の認識を、一chosen pointから全源へ根拠なく外挿せず、ここで証明したwhole closure等式で全域へ戻せる。

現時点では実軌道次元は旧上界23以下だけを使う。新ambient D、成分性、親、相互相違は本純証明だけでは未判定であり、後続の精確な有限義務とする。
