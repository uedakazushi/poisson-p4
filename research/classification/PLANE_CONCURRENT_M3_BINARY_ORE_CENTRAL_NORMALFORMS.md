# 全binary Ore Qのcentral二次項：全rank壁とactual正常化

2026-11-08、PRIMARY G、producer純証明・独立review前。
新CAS、source probe、minor farm、GB又はjetは無い。
先行全Q coverのC2三係数を可逆なactual changesで正常化する。
全unclosed familyはmagneticを含めて高々二parametersとなるが、full W moduleのG上界は未了。

## 1. 共通shift mapと真のmagnetic scale

全Qは

`Qab=a A1+b B1+C2`, `Qan=U`, `Qbn=V`, `Quv=0`。

actual a'=a+ell_a(u,v),b'=b+ell_b(u,v)は

`C2'=C2+U(ell_b)−V(ell_a)−A1 ell_a−B1 ell_b`。

ell_a=u,v、ell_b=u,vの順に、全shift imageは四columns

```
S1=−V(u)−A1u, S2=−V(v)−A1v,
S3= U(u)−B1u, S4= U(v)−B1v.                          (1)
```

に尽くされる。全inverseはminus shifts。
normal variablesを固定しa,bをλ^-1倍、full tensorをλ倍すれば、
U,V,A1,B1は元のliteral値に戻り、C2だけをλ^-1倍する。
従って残一coefficientはzero/nonzeroに分け、非零を1へactualに戻せる。
このscaleはfull Wも変えるが、全quadratic W sourceへ戻る可逆変換である。
残二coefficientsを独立scale出来るとは自動的に仮定しない。

## 2. K非零の全magnetic rank strata

先行noteのcanonical symbols/parametersをそのまま使う。

### K211

S1=S3=S4零、S2=3v(u−v)。全tでrank一。
imageを除くとC2=x u²+y v²である。
共通magnetic scaleにより、全supportはC2零、又は全projective P1の二charts

`u²+z v²`, `z u²+v²`

で覆う。tとmagnetic projective coordinateを保持するので最大二parameters。
root marksを保つ任意でないnormal torusを使いx,yを同時に1へはしない。

### K22

全columnsは

`S1=S4=0, S2=(β1−β2)v², S3=(α1−α2)u²`。

det(α,β)≠0だからα1=α2とβ1=β2が同時には起きない。
両差非零ではimage span(u²,v²)、residue C2=γuv、γ0/1。
α1=α2ではimage span(v²)、residues(u²,uv)。
β1=β2ではimage span(u²)、residues(uv,v²)。

後二場合のresidual二coefficientsも真に独立scale出来る。
normal u'=s u,v'=t v、base a'=s a/λ,b'=t b/λ、tensor scaleλはcanonical α,βを保ち、
C2の(u²,uv,v²)coefficientsを

`(t/(λs), 1/λ, s/(λt))`

倍する。各residual pairの二charactersは独立なので各zero/nonzeroの全四flagsへ移る。
全P1 rows/chart boundariesとdet≠0を保持する。
genericは二parameters、rank一wallは一parameter以下。

### K31

C=[[r,s],[t,−r]],r²+st=0について

```
S1=(t−r)uv, S2=−uv+2t v²,
S3=(s−r)uv, S4=−uv+(2s−2r)v².
```

C非零ではrank二、image=span(uv,v²)。
実際t−r=s−r=0ならnil conditionが全r,s,t零を強制する。
少なくとも一つからuvを得、v²が得られない場合も同じnil conditionでC零となる。
residueはγu²、γ0/1。
C零ではimage span(uv)、residues u²,v²。
このsource U=(u²,0),V=(0,uv)はnormal u'=s u,v'=t v、base a'=s a/λ,b'=s b/λ、tensorλで保たれ、
u²,v²のcoefficientsは(1/λ,s²/(λt²))倍となる。
複素非零scalesで全四0/1 flagsへ移る。CのP1/zero全coverを保持する。

### K4a

全columnsは

```
S1=−αu²+2(1−β)uv,
S2=−u²−2αuv+3(1−β)v²,
S3=(β−1)u², S4=αu²+2(β−1)uv.
```

β≠1ではrank三でC2零へ戻せる。
β1,α≠0ではrank二image span(u²,uv)、C2=γv²、γ0/1。
β1,α0ではrank一image span(u²)、residuals uv,v²。
このspecial pairはU=uE_n,V=vE_n+u²∂v。
actual v'=v+k u、a'=a,b'=b+k aはpair/A1/B1を保ち、
residual C2=xuv+yv²を `(x−2ky)uv+yv² mod u²` へ動かす。
y非零ならxを除け、共通mag scaleでy1。
y零ならx0/1。従って全補集合は `0,uv,v²` の三flags。

### K4b

全columnsは

```
S1=0, S2=−u²,
S3=(α−1)u²+2βuv, S4=2(α−1)uv+3βv².
```

β非零rank三、C2零。
β0,α≠1はrank二、residueγv²、γ0/1。
β0,α1はrank一image span(u²)。
このpairU=uE_n,V=u²∂vはv→v+k uだけで保たれる。
同じresidual shearとcommon scaleにより全三flags0,uv,v²。
全rank一/二の補集合も元の二parameter族へactualに戻る。

## 3. K零constant-direction全四類

### C_tan

```
S1=−u²−tuv, S2=−2uv−t v²,
S3=t u², S4=u²+tuv.
```

t非零rank三、C2零。t0rank二image span(u²,uv)、residueγv²、γ0/1。

### C_trans

imageはdiagonalに

`(t+3)u², (t+2)uv, (t+1)v²`

である。
tが−3,−2,−1以外ならrank三。
各三wallはrank二でmissing monomialは順にu²,uv,v²。
各residueγ0/1を保持する。

### B_axis/B_other

B_axisは(S1,S2,S3,S4)=(-v²,0,3u²,2uv)。
B_otherは(-v²−2uv,-3v²,3u²,u²+2uv)。
どちらもconstant rank三、全C2零へ戻せる。

## 4. K零SS/Jordanの全period壁とcenter空性

SS(t,s)、t≠0のimageは

`(s+2)u², (s+t+1)uv, (s+2t)v²`。

三affine linesは一点(t1,s−2)でconcurrentである。
三lines外はrank三、C2零。
各lineからcenterを除けばrank二で、そのlineのmissing monomialγ0/1を保持する。
各line全parameterとt0の別C_trans source complementを保つ。

centerではshift rank零だが、C2三係数を捨ててはいけない。
ここでは

```
U=uE_n, V=vE_n,
A1=−v, B1=u,
Q=−E_y∧(u∂a+v∂b)+C2∂a∧∂b.
```

PfQ零、actualD=−[E_y∧Q∧W] はC2 times homogeneous cubic。
C2非零ならproper factor、C2零ならD零。
したがって**centerの全三magnetic sourceは空**、Wのchosen pointでのfactorではない。

Jordan(s)のcolumnsは

```
S1=−u²−(s+2)uv, S2=−2uv−(s+2)v²,
S3=(s+2)u², S4=u²+(s+2)uv.
```

s≠−2でrank三、C2零。s−2でrank二、residueγv²、γ0/1。
nil/zero-eigenvalueのC_tan補集合は先行whole Q proofの通りである。

## 5. 全sourceの残parameter数と未了義務

全K≠0 sourceは

* K211：t＋全projective magnetic P1、又はzero magnetic。
* K22：ordered projective rowsの全二parameter generic0/1、及び全一parameter row-equality wallsの四flags。
* K31：projective nilcone P1の0/1 flags、C零の四flags。
* K4a/b：全二parameter C2零、全一parameter wallsの0/1、special pointの三flags。

全K0 sourceは

* C_tan/C_trans：全tのC2零、有限t wallsの0/1。
* B_axis/B_other：fixed C2零。
* SS：全二parameterC2零、三全一parameter wallsの0/1、centerは全空。
* Jordan：全sのC2零、s−2の0/1。
* constant-proportional U,V：既受理global essential三上界23。
* U=V零：全actualD factor又は零。

で尽くされる。
全unclosed canonical Qは**centralを含め高々二parameters**。
全mapsはactual linear block/scalar/positive shearで、zero/nonzeroを明示し逆を保持した。
full W moduleは依然全quadratic WとE_y∧L_WQ0、元normalweight/actualD3/R conditionsである。
この全moduleの像次元、parameter rank walls、normalizer fiber又はactual factorが次の必要最弱義務。
今回の小さいshift mapの正常化をfullPoisson分類完了と呼ばない。
