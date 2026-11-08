# R=3L・m二SF02全空性の独立手証明：係数三とcurl係数二

2026-11-07、PRIMARY G、`/root/log_manuscript`。
`PLANE_TRIPLE_LINE_M2_SQUAREFREE_FLAG_02_EXCLUSION.md`
（SHA256 `e9b6c1ad23d413ee80333ff31bf929d6d7a53b1b18bd02c008adcdbc5f336cc7`）を
別作者として全読解し、全係数・sign・complementsを以下の手導出で確認した。**PASS**。
新point、probe、matrix、CAS、jet又はGBは実行していない。
新成分・他flags・全G/Eの完了は主張しない。

## 1. 全irregular connectionとtracefree fields

actualD被約既約quintic、rank二restriction plane、f=a³、normal位数二、generic SF、flag02を仮定する。
`PLANE_NONREDUCED_M2_SPLIT_ROOT_ENTRANCE.md` の全pole budgetを再読し、
primitive rowのunitからcharacter pole上界を先に得ること、
`pole_divisor(eta)≤div(h)` が非被約collisionも含むことを確認した。

全root sourceは非零amplitude gammaを保持して

```
G2=gamma*a*u*(A*u+a²v), gamma!=0,
M=[[1,0],[A,a²]], gcd(A,a)=1,
eta=−2dlog a.
```

直接M逆を戻すと `beta11=beta22=theta`、beta12零、
`beta21=−d(A/a²)` である。
thetaは `dlog a+d(P2/a²)`、全P2 modulo Ca²の五exact係数を含む。
これをlog部分だけに制限していない。

全P0はP_f、mixed columns `uX` と `uZ+vX`、normal bracket
`k u²+kappa uv+k02 v²` の形で、X,Zはtracefree linear。
intrinsic c一chartのanchor3fとEuler補正を戻すと

```
alpha_X=iE iXvol/f=dlog f−3theta,
iE iZvol=−3f beta21=3a dA−6A da.
```

closed thetaとtrace X零に対する直接Cartan式
`d alpha_X=X(f)iEvol/f²` からX(f)=0、従ってXa=0。
第二式は `Z=(0,−3A_c,3A_b)` を与える。
実際このZに対するiE iZvolは
`3a(A_b db+A_c dc)−3(bA_b+cA_c)da=3a dA−6A da`。
同じcontractionのkernelはEuler multipleで、tracefree条件によりそれも零である。
したがって全a-column translationsと全trace零2×2 block Hを保持する。

## 2. actual normal coefficientsのsignとfactor三

座標順(a,b,c,u,v)でactualDを `E_total∧P∧P/2` とする。
任意のX,Y,Zに対するnormal二次係数は

```
u²: 3f k−det(E,X,Z)−Z(f),
uv: 3f kappa−det(E,X,Y)+X(f)−Y(f),
v²: 3f k02.
```

第一の3はEulerのcubic identity E(f)=3f、
−Z(f)はEulerのu項、+X(f)はEulerのv項から直接生じる。
今はY=X、Za=0なのでcross termsは零、v²比較でk02零、uv比較で
**gamma=3kappa、kappa非零**。
さらに上のZから `det(E,X,Z)=3a X(A)` だから、u²比較は

```
X(A)=a²k−kappa A.                                (1)
```

全P0のactual Pfaffian consistencyであり、flat connectionが存在するだけでP0を仮定する議論ではない。

## 3. 元projective一式からcurl係数二対三

全raw P1三十六/Euler六と全P2九を保持する。
u-columnのbad部分は `v²C`、C constant。
P2にnormal columnは無く、C零なら全u-rowがu倍となってactualDにu因子を与える。
従って指定全sourceではC非零が必要である。

u零でstrict cyclic Jacobiの係数を直接展開すると

```
J_iuv[v³]=kappa C_i+(DY C)_i,
J_iju[v²]=[(2Y−X) wedge C]_ij−(P_{C(f)})_ij.
```

前者のkappa項は `P_iu*partial_u P_uv`、DY項は
`P_uj*partial_j P_vi` で正符号。
後者のmixed-v derivativeが2Y、mixed-u derivativeが−X、base derivativeが−P_Cfである。
Z/kはu因子を持ち、他P1/P2は別normal monomialなので寄与しない。
strict Jacobi自体を零とは置かず、元quad(i,j,u,v)のprojective条件だけから

```
E wedge(kappa C+DY C)−(2Y−X) wedge C+P_{C(f)}=0.
```

base volume双対でcurlを取れば
`curl(E×D)=−2D`、`curl(V×C)=DV C`（trace V零）、gradient項のcurl零なので

```
(2kappa I+4DY−DX)C=0.
```

この導出はR被約性を一度も使わず、f=a³にもそのまま適用できる。
Y=Xの本源では正確に

```
(2kappa I+3DX)C=0.                               (2)
```

## 4. 全translationsとnonsemisimple境界を残す矛盾

Xa=0なのでDXのa-rowは零で、(2)からCa零。
C非零だから残りの二成分はHの固有vector、固有値は−2kappa/3。
trace H零とkappa非零からHの二固有値はlambda,−lambdaで
`lambda=±2kappa/3`、相異なる非零値である。
nilpotent又はzero Hはここで除外済みであり、初めからgeneric Hを仮定していない。

(1)をa零へ制限すると、primitivityにより非零binary quadratic
`Abar(b,c)=A(0,b,c)` が
`H(Abar)=−kappa Abar` を満たす。
Sym²の固有値は2lambda,0,−2lambdaなので
`lambda=±kappa/2` が必要。
±2/3と±1/2は複素数上も一致しないから矛盾する。
この証明はa-columnの二translationsを零と置かずに成立する。
望むならH可逆性でactual shearにより消せるが、判定には不要である。

以上でC非零が不可能、C零も既約actualDに反するため、**全flag02 sourceは空**。
全theta exact coefficients、primitive Aの全係数、codimension二basepoints、
全P1/P2とnonzero amplitudeを保持した全点空性である。
f=a²b、flags00/01、normal double quadratic、m三又は全nonreduced-m二への拡張はしていない。
