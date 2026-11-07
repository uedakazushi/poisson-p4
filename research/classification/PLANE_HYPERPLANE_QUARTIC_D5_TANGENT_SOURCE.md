# D5 quartic と cusp tangent：全marked family・全六磁気源の一軌道還元

2026-11-07、PRIMARY E、producer `exceptional_coverage`。手による全族証明、独立確認前。
対象は integral quartic Fが唯一D5を持ち、Hがそのcuspidal branchの接線、H∩Fはその一点だけ、
かつD=HFがexponents(2,2)のfree quinticとなるmarked branchである。
他のD5 marked line、他の九profile、全Hquartic又は全Eの尽くしは主張しない。

## 1. 全marked quarticを一つの正規形へ

H=(h=0)、D5点p=[0:0:1]、transverse smooth branchの接線をy=0に取る。
全quarticは、非零scalarを除き

```
F=h²*y*w+a*h⁴+b*h³*y+c*h²*y²+d*h*y³+e*y⁴
```

である。multiplicity三とtangent cone h²yを全て使っただけで、残るquartic五係数を保持している。
Hが成分でないのでe!=0。可逆なw-shearでb,cを消せる。
従って全marked familyをF=h²yw+a h4+d h y3+e y4とする。

affine w=1のg=hFは、weights(h,y)=(3,2)に対して

```
g0=h³*y+e*h*y⁴, wt(g0)=11,
g=g0+d*h²*y³+a*h⁵,
wt(h²*y³)=12, wt(h⁵)=15.
```

g0の二partialはleading weights八、九のregular sequenceである。
そのweighted Milnor algebraのHilbert seriesは
`(1-t^8)(1-t^9)/((1-t^3)(1-t^2))`、length十二、top weight十二。
gのpartialのinitial formsは同じregular sequenceであるから、weighted filtrationで
gr(Q_g)=Q_g0。これはinitial syzygyがKoszulだけであり、逐次liftできることから従う。
局所formal completionで証明し、有限lengthのMilnor/Tjurina代数ではanalyticと一致する。

weight十二の単項式は `(y⁶,h²y³,h⁴)`。
Jacobian idealの同weight部分は正確に

```
y²*(g0)_h=3h²y³+e*y⁶,
h*(g0)_y=h⁴+4e*h²y³
```

の二次元包であり、h²y³のclassは非零。
weighted Euler identityをQ_g内で使うと

```
-11[g]=d*[h²y³]+4a*[h⁵].
```

h5はtop weight十二より高いのでQ_g内で零。
d!=0なら右辺のinitial weight十二は非零で、gはJacobian idealに入らずtau(g)<mu(g)=12。
d=0ならgはJacobian idealに入りtau(g)=mu(g)=12。
FのD5はdelta三を使い切り、Hとの交点もpだけなので、これは全Dのtau条件そのもの。
従って全対象で **d=0** が必要十分である。

d=0でa=0ならF=y(h²w+e*y³)が可約なのでa!=0。
a,e!=0ならwの一次式としての二係数h²yとa h4+e y4は互いに素で、Fは絶対既約。
対角base変換とquartic scalarでa=e=1に直る。具体的にh=H,y=beta*Y,w=(a/beta)*W、
beta⁴=a/eとすればF=a*(H²YW+H⁴+Y⁴)。全行列は可逆で、逆も保持できる。
以上で全marked familyは一つのpair

```
H=h, F=h²*y*w+h⁴+y⁴
```

のPGL3軌道に入る。Q7という一代表だけから全familyを推測した議論ではない。

## 2. 全二次Saito anchorと実Poisson源

h=1でF=yw+1+y4とする。全quadratic polynomial vector field
U=(U1,U2)についてU(F)=ell*F、ellは全affine linear、と置いて係数比較すると、全解は

```
U = A*(y²,4+3yw)+B*(4+yw,-16y²+3w²).
```

ここでA,Bは定数である。全二成分各六係数とell全三係数を比較した全空間の解であり、
一つのsyzygy ansatzではない。基底を

```
U=y²*partial_y+(4+3yw)*partial_w,
V=(4+yw)*partial_y+(-16y²+3w²)*partial_w
```

とすると、直接

```
U(F)=4y*F, V(F)=4w*F,
det(U,V)=-16F,
[U,V]=3w*U-y*V.
```

homogenizeしたE,U,Vのdeterminantは−16hF。よって元のfree2,2条件を実際に満たす。
任意の全二次magnetic form Rについて

```
P=partial_q1 wedge U+partial_q2 wedge V
  +(3w*q1-y*q2+R(y,w))*partial_q1 wedge partial_q2
```

は全Jacobiを満たす。上のcommutator identityを代入すれば直接分かる。
Pfaffianは16F、実五次は16hF。全sourceの六磁気項を一つも捨てない。
既済の全cotangent/Saito入口により、任意のsourceframeは定数GL2 fiber変換でこの基底へ移り、
base変更の二次形式pullbackも可逆なので、元の全marked family・全六磁気を含む。

## 3. 全六磁気の多項式正規化と零境界

R=m0+m1*y+m2*w+m3*y²+m4*yw+m5*w²とする。
実可逆変換q1'=q1+L1,q2'=q2+L2について
R'=R−3wL1+yL2+U(L2)−V(L1)である。
全六affine L係数の作用は、順に

```
L1=1,y,w: -3w, -4-4yw, 16y²-6w²,
L2=1,y,w: y, 2y², 4+4yw.
```

rank五、唯一のcokernel invariantは I=m0−m4。
全parameterで有効な具体的normalizerは

```
L1=m2/3+(m5/6)w,
L2=-m1-(m3/2+4m5/3)y-(m4/4)w,
R'=I.
```

磁気係数は一つも反転しない。I!=0ならq_new=q/I,tensor multiply IでI=1へ移る。
I=0の全補集合は、正則な実族R=tのt=0端点である。
従って全六mag sourceのPGL飽和閉包は、canonical R=1の一つの実軌道の閉包に等しい。

## 4. canonical軌道は正確に23次元

marked base pairのinfinitesimal PGL3 automorphismは零である。
唯一のD5点とcuspidal/smooth二接線を保つので、linear fieldは
X_h=alpha*h,X_y=gamma*y,X_w=delta*h+epsilon*y+zeta*w。
X(F)=kappa*Fを比較するとdelta=epsilon=0、alpha=gamma=zetaとなり、projectively scalarだけ。

実D=16hFのessential rankは三で、そのdirectional-derivative kernelはfiber二平面である。
全PGL5 stabilizerはこのintrinsic kernelと唯一の一次因子Hを保ち、base quotient上は上のfinite groupへ写る。
従ってinfinitesimallyはq'=Aq+L(h,y,w)、base恒等だけを調べればよい。
projective tensor equality G_*P=sigma*Pのanchor成分はA*M=sigma*MなのでA=sigma*I。
磁気cokernel I=1はsigma²=sigmaを強制する。sigma!=0なのでsigma=1。
残るtranslation kernelは正確にL1=y,L2=wの直線、field N=y*partial_q1+w*partial_q2だけ。
ゆえにprojective stabilizer dimensionは一、PGL5軌道dimensionは24−1=23。
全source閉包が一つの23次元軌道閉包であることを、fullambient行列計算前に証明した。

## 5. まだ必要な検証

この23次元source閉包が全Poisson成分か、より大きい既知成分の境界かは未解決。
正当化された唯一のcanonical R=1点で全150/350 tangentと真正cutを調べる有限jobを次に行う。
全D102ならsmooth23、D101なら一本のnormalの最弱非零二次余核で十分。
D<=100なら深いjetへ進まず、その弱い結論と障害を保存する。
既存全成分との相違も、成分性の後に閉因子・essentialrank・marked配置で全件監査する。
本手証明だけではcomponent数を増やさない。
