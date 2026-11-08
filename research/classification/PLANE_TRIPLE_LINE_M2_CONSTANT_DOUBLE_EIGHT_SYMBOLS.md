# R=3L・m二constant double：全源は八つの固定P0又は既約性矛盾

2026-11-08、root純手証明。独立review前、higher算術未実行。
対象はactualD被約既約quintic、rank二restriction plane、f=a³、normal位数二、
generic double normal quadraticでprimitive root degree零の全枝である。
全irregular/off-diagonal接続を含む必要還元を示す。全G/E又は八つのhigher源の閉鎖はまだ結論しない。

## 1. full P0と元projective bad-vector式

constant normal右GL2でG2=gamma*a³*u²、gamma非零とする。
parallel leading relationはdiagonal ordinary connectionの差を零にし、
全closed diagonal、off-diagonalを戻すとtracefree mixed gaugeは

```
Pbase=P_(a³), {base,u}=uX, {base,v}=uZ+vX,
{u,v}=k*u².
```

closedness/flatnessからX(f)=Z(f)=0、従ってXa=Za=0。
全irregular potentialsはX,Zの全tracefree linear coefficientsへ含まれ、先にlog部分へ限定していない。
aがP0のCasimirなのでprojective Jacobiはstrict Jacobiと同値、全式から[X,Z]=0。
tracefree二fieldsはa=1上でaffine Hamiltonian fieldsで、det(X_bc,Z_bc)=d0*a²はconstant Hamiltonian bracket。
actualG2は(3k−d0)*a³*u²、gamma=3k−d0非零である。

全raw36/Euler六・全P2九を許す。
u-columnのbad v²CはEuler gaugeで変わらず、C零ならfullactualDにu因子が生じるためC非零。
元projective weight一の全bad-vector恒等式は

```
E wedge(DX*C)−X wedge C+P_(C(f))=0.
```

そのcurlは3DX*C零を与える。従ってDX*C零で、上式は
`−X wedge C+P_(3Ca*a²)=0`。
X非零ならab/ac成分がCa*Xb=Ca*Xc零を要求し、Ca零。
続いてX wedge C零となり、Xの全imageは一つのconstant vector Cのspanに含まれる。
従ってX非零ならrank DX一かつ(DX)²零、Xa零である。
semisimple X又はregular nilpotent三block Xは、この全higher必要式だけで空となる。

## 2. k非零・nilpotent Zは全higherで空

これは既 `PLANE_M2_DOUBLE_EQUAL_DIAGONAL_NIL_OFF_EXCLUSION.md` の証明を、
leading係数3kの代わりに任意gamma非零へそのまま適用する場合である。
仮定としてk非零、X=Y、Za=0、Z nilpotent、actualF2=gamma*a³*u²を用いる。
det(E,X,Z)零はこの段階に不要であり、k非零をgammaから推論してもいない。

canonical strict Casimir liftのM=Pi sharp d(v/u)は
Mnormals=3k(u,v)/5、Mbase=(Z−2kI/5)base+B(u,v)となる。
kI−Zはunit determinantなので、normal degree一の有理translationでBを消せる。
base degree5−mの全summandへのoperatorはZ+k(m−2)I。
nilpotenceによりm二以外可逆、従ってexactly

`F=gamma*u²*(a+T_a(u,v))³`。

C(u,v)上でpositive base-degree factorsを持ち、Gauss lemmaからpolynomial Fの既約性に反する。
全P1三十/P2九を許した証明である。有理translationをglobal PGLとは呼ばない。

## 3. X非零の全branches：二つのnil固定源だけ

rank一nilpotent Xはbase GL3でa-lineを保ちT=partial_cとして

`X=(r*a+s*b)*partial_c`

と書ける。r,sは同時零でない。
全Za零/trace Z零のlinear Zをcommuting式へ戻す。
s非零ならそのb-componentが零と強制され、Z=(C*b+q*a)*partial_cである。
従ってd0零、gamma=3k非零、Z nilpotentなので前節で全空。

s零ならr非零、全commuting Zは

`Z=p*a*partial_b+(C*b+q*a)*partial_c`。

これもnilpotent。k非零は前節で空、k零ならgamma=r*p非零なのでp非零。
base c→c+ell*bの可逆shearでqを消せる。
anchor3a²を保つbase diagonal/tensor scaleとrelative normal u/v scaleで、
r,pを一、Cを零又は一に合わせられる。具体的にはold b=beta*b,old c=chi*c,
tensor sigma=beta*chiならr'=beta*r、p'=e*chi*p、C'=e*beta²*C。
beta=1/rを固定し、C非零はe=1/(beta²*C)、chi=beta²*C/p、
C零はe一、chi=1/pと取る。全unitとその逆を保持する。

従って残る全sourceは

```
X=Y=a*partial_c,
Z=a*partial_b+epsilon*b*partial_c,
k=0, epsilon=0,1, gamma=1
```

の二固定P0へ入る。normal shearで勝手にparameterを消したのではない。

## 4. X零の全branches：六つの半単純固定源だけ

X零ならgamma=3k非零、k非零。全tracefree Za零のlinear Zを保つ。
bc blockの固有値が零ならfull Zはnilpotentで、§2により空。
相異なる固有値lambda,-lambda（lambda非零）ならactual base conjugationとa-column shearで

`Z=lambda*(b*partial_b−c*partial_c)`

へ戻せる。相対normal scale/tensor scaleでk一にし、lambda/kを新lambdaとする。
全higher canonical ratio式でA=Z、translationはkI−Zが可逆なら消せる。
normal degree m=3/4/5、base degree d=5−mの必要operatorは
Sym^d(Z)+k(m−2)Iである。
base spectrumはn*lambda、n=-d,-d+1,...,d。
従ってhigher correctionを許す全spectral wallsは

```
m=3: lambda/k=±1,±1/2,
m=4: lambda/k=±2,
m=5: none.
```

translation failure lambda/k=±1もこの表へ含まれる。
これら以外ではF=gamma*u²*(a+T_a)³で既約性に反する。
lambda零/nilpotent及びk零は既に別branchで処理し、全分母壁を捨てていない。
残る全X零源は

```
X=Y=0, Z=lambda*(b*partial_b−c*partial_c), k=1,
lambda∈{1/2,-1/2,1,-1,2,-2}, gamma=3
```

の六固定P0へ入る。

## 5. 次の義務

二nil固定＋六sem固定の八つで全first-symbol源を必要被覆した。
元のfullhigher coefficient三十/九は残っており、次の最弱dimension/factor検査なら八whole sourcesを一回で扱える。
本noteはそのG上界又はall-point E帰属を証明したとは呼ばない。
R=a²b、m二SF、moving double、m三等の別源もこの八表から尽くしたとは主張しない。
全connectionからXa=Za=0及び[X,Z]=0、bad-vector式のfactorとspectral壁のscopeは独立reviewを要する。
